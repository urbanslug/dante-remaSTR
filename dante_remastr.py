from datetime import datetime
from typing import TextIO, Iterable, Iterator, Any
import argparse
import csv
import gzip
import io
import itertools
import multiprocessing
import re
import os
import sys

import numpy as np
import pandas as pd
import numpy.typing as npt

from src.inference import Inference, phase, save_phasing
from src.arguments import load_arguments
from src.report import (
    chrom_from_string, ChromEnum, write_all, log_str, write_alignment,
    write_histogram_nomenclature, configure_logger, write_report
)
from src.annotation import (
    Annotation, pairs_to_annotations_pick, annotations_to_pairs,
    remove_pcr_duplicates
)
from src.annotation.Motif import Motif
from src.postfilter import PostFilter
from src.filtering import has_good_quality, cut_low_quality


def shorten_str(string: str, max_length: int = 40, ellipsis_str: str = '...') -> str:
    """
    Shorten string to max_length and include ellipsis.
    :param string: str - string to shorten
    :param max_length: int - maximum length to shorten to
    :param ellipsis_str: str - ellipsis to add to end of string
    :return: str - shortened string
    """
    if len(string) > max_length:
        return string[:max_length - len(ellipsis_str)] + ellipsis_str
    else:
        return string


def errors_per_read(
    errors: list[tuple[int, int, int]], relative: bool = False
) -> tuple[float | str, float | str]:
    """
    Count number of errors per read. Relative per length or absolute number.
    :param errors: list[tuple[int, int, int]] - indels, mismatches and length of module
    :param relative: bool - relative?
    :return: tuple[float, float] - number of indels, mismatches per hundred reads
    """
    # if we have no reads, return '---'
    if len(errors) == 0:
        return '---', '---'

    if relative:
        mean_length = np.mean([length for _, _, length in errors])
        return (
            float(np.mean([indels / float(length) for indels, _, length in errors]) * mean_length),
            float(np.mean([mismatches / float(length) for _, mismatches, length in errors]) * mean_length)
        )
    else:
        return (
            float(np.mean([indels for indels, _, _ in errors])),
            float(np.mean([mismatches for _, mismatches, _ in errors]))
        )


def generate_result_line(
    motif_class: Motif, predicted: tuple[str | int, str | int], confidence: tuple[float | str, ...],
    qual_num: int, primer_num: int, filt_num: int, module_number: int,
    qual_annot: list[Annotation] | None = None,
    second_module_number: int | None = None
) -> dict:
    """
    Generate result line from the template string.
    :param motif_class: Motif - motif class
    :param predicted: tuple[str, str] - predicted alleles (number or 'B'/'E')
    :param confidence: tuple[7x float/str] - confidences of prediction
    :param qual_num: int - number of reads with both primers
    :param primer_num: int - number of reads with exactly one primer
    :param filt_num: int - number of filtered out reads (no primers, many errors, ...)
    :param module_number: int - module number in motif
    :param qual_annot: list[Annotation] - list of quality annotations for error and number of reads
    :param second_module_number: int/None - second module number in motif
    :return: dict - result dictionary
    """
    # setup motif info
    start, end = motif_class.get_location_subpart(module_number)
    motif_seq = motif_class.module_str(module_number)
    if second_module_number is not None:
        _, end = motif_class.get_location_subpart(second_module_number)
        motif_seq = ','.join(
            [motif_class.module_str(i) for i in range(module_number, second_module_number + 1)]
        )

    reads_a1: int | str
    reads_a2: int | str
    indels_rel: float | str
    indels_rel1: float | str
    indels_rel2: float | str
    mismatches_rel: float | str
    mismatches_rel1: float | str
    mismatches_rel2: float | str
    # get info about errors and number of reads from quality annotations if provided
    reads_a1 = reads_a2 = '---'
    indels_rel = mismatches_rel = '---'
    indels_rel1 = mismatches_rel1 = '---'
    indels_rel2 = mismatches_rel2 = '---'
    if qual_annot is not None:
        # get info about number of reads
        a1 = int(predicted[0]) if isinstance(predicted[0], int) else None
        a2 = int(predicted[1]) if isinstance(predicted[1], int) else None
        reads_a1 = 0 if a1 is None else len(
            [a for a in qual_annot if a.module_repetitions[module_number] == a1]
        )
        reads_a2 = 0 if a2 is None else len(
            [a for a in qual_annot if a.module_repetitions[module_number] == a2]
        )

        # get info about errors
        errors = [a.get_module_errors(module_number) for a in qual_annot]
        errors_a1 = [a.get_module_errors(module_number) for a in qual_annot
                     if a.module_repetitions[module_number] == a1]
        errors_a2 = [a.get_module_errors(module_number) for a in qual_annot
                     if a.module_repetitions[module_number] == a2]
        assert len([l for i, m, l in errors if l == 0]) == 0

        # extract error metrics
        indels_rel, mismatches_rel = errors_per_read(errors, relative=True)
        indels_rel1, mismatches_rel1 = errors_per_read(errors_a1, relative=True)
        indels_rel2, mismatches_rel2 = errors_per_read(errors_a2, relative=True)

    # return dictionary
    return {
        'motif_name': motif_class.name, 'motif_nomenclature': motif_class.motif, 'motif_sequence': motif_seq,
        'chromosome': motif_class.chrom, 'start': start, 'end': end, 'allele1': predicted[0],
        'allele2': predicted[1], 'confidence': confidence[0], 'conf_allele1': confidence[1],
        'conf_allele2': confidence[2], 'reads_a1': reads_a1, 'reads_a2': reads_a2, 'indels': indels_rel,
        'mismatches': mismatches_rel, 'indels_a1': indels_rel1, 'indels_a2': indels_rel2,
        'mismatches_a1': mismatches_rel1, 'mismatches_a2': mismatches_rel2, 'quality_reads': qual_num,
        'one_primer_reads': primer_num, 'filtered_reads': filt_num,
        'conf_background': confidence[3] if len(confidence) > 3 else '---',
        'conf_background_all': confidence[4] if len(confidence) > 4 else '---',
        'conf_extended': confidence[5] if len(confidence) > 5 else '---',
        'conf_extended_all': confidence[6] if len(confidence) > 6 else '---',
        'repetition_index': module_number
        if second_module_number is None
        else f'{module_number}_{second_module_number}'
    }


# this has mixed functionality - it does prediction and it does writing
def too_complex_function_inferring_and_creating_reports(
    args: argparse.Namespace,
    annotations1: list[Annotation],
    postfilter_class: PostFilter,
    motif_class: Motif,
    motif_dir: str,
    motif_str: str,
    module_number: int,
    read_distribution: npt.NDArray[np.int64],
    result_lines: list,
    prev_module: tuple[int, str, int] | None,
    file_pcolor: str | None,
    file_output: str | None
):
    monoallelic1 = args.male and chrom_from_string(motif_class.chrom) in [
        ChromEnum.X, ChromEnum.Y
    ]

    # setup post filtering - no primers, insufficient quality, ...
    qual_annot1, filt_annot1 = postfilter_class.get_filtered(annotations1, module_number, both_primers=True)
    primer_annot1, filt_annot2 = postfilter_class.get_filtered(filt_annot1, module_number, both_primers=False)

    # why is this a class?
    # run inference - this takes most of the time (for no --verbose)
    model = Inference(
        read_distribution, args.param_file, str_rep=args.min_rep_cnt,
        minl_primer1=args.min_flank_len, minl_primer2=args.min_flank_len,
        minl_str=args.min_rep_len
    )

    predicted1, confidence1 = model.genotype(
        qual_annot1, primer_annot1, module_number, file_pcolor, file_output,
        motif_str, monoallelic1
    )

    # append to the result line
    result_lines.append(generate_result_line(
        motif_class, predicted1, confidence1,
        len(qual_annot1), len(primer_annot1), len(filt_annot1),
        module_number, qual_annot=qual_annot1
    ))

    # infer phasing (if we are not on the first repeating module)
    last_num1 = None
    both_good_annot1 = None
    one_good_annot1 = None
    none_good_annot1 = None
    phasing1 = None
    supp_reads1 = None
    if prev_module is not None:
        # get the last module number
        last_num1 = prev_module[0]

        # post filtering
        both_good_annot1, filtered_annot1 = postfilter_class.get_filtered_list(
            annotations1, [last_num1, module_number], both_primers=[True, True]
        )
        left_good_annot1, left_bad_annot1 = postfilter_class.get_filtered_list(
            filtered_annot1, [last_num1, module_number], both_primers=[False, True]
        )
        right_good_annot1, none_good_annot1 = postfilter_class.get_filtered_list(
            left_bad_annot1, [last_num1, module_number], both_primers=[True, False]
        )
        one_good_annot1 = left_good_annot1 + right_good_annot1

        # infer phasing
        phasing1, supp_reads1 = phase(both_good_annot1, last_num1, module_number)

        # append to the result line
        result_lines.append(generate_result_line(
            motif_class, phasing1, supp_reads1,
            len(both_good_annot1), len(one_good_annot1), len(none_good_annot1),
            last_num1, second_module_number=module_number
        ))

    return (
        (qual_annot1, primer_annot1, filt_annot2),
        (both_good_annot1, one_good_annot1, none_good_annot1),
        phasing1, supp_reads1, confidence1, predicted1, result_lines, last_num1
    )


def report(
    args, motif_dir, motif_class, module_number, last_num1, prev_module,
    ann_qual, ann_primer, ann_filter, ann_2good, ann_1good, ann_0good,
    phasing, supp_reads, confidence, predicted
):
    # TODO: merge meaning prev_module and last_num1

    # write files if needed
    write_all(
        ann_qual, ann_primer, ann_filter,
        motif_dir, motif_class, module_number,
        zip_it=args.gzip_outputs, cutoff_alignments=args.cutoff_alignments
    )

    if prev_module is not None:
        # write files
        write_all(
            ann_2good, ann_1good, ann_0good,
            motif_dir, motif_class, last_num1, module_number,
            zip_it=args.gzip_outputs, cutoff_alignments=args.cutoff_alignments
        )

        # write phasing into a file
        save_phasing(
            f'{motif_dir}/phasing_{last_num1}_{module_number}.txt', phasing, supp_reads
        )

    # write the alignments
    if confidence is not None:
        # get number of precise alignments for each allele
        a1 = int(predicted[0]) if isinstance(predicted[0], int) else None
        a2 = int(predicted[1]) if isinstance(predicted[1], int) else None

        if a1 is not None and a1 > 0:
            write_alignment(
                f'{motif_dir}/alignment_{module_number}_a{a1}.fasta',
                ann_qual, module_number, a1,
                zip_it=args.gzip_outputs, cutoff_after=args.cutoff_alignments
            )
        if a2 is not None and a2 != a1 and a2 != 0:
            write_alignment(
                f'{motif_dir}/alignment_{module_number}_a{a2}.fasta',
                ann_qual, module_number, a2,
                zip_it=args.gzip_outputs, cutoff_after=args.cutoff_alignments
            )


def process_group(
    args: argparse.Namespace, df: pd.DataFrame, motif_str: str
) -> tuple[Motif, list[dict], int, int]:
    """
    Process the group as pandas Dataframe. Return motif name if processed correctly or None otherwise.
    :param args: argparse.Namespace - namespace of program arguments
    :param df: pd.DataFrame - contains information about annotated reads for a single motif to process
    :param motif_str: str - motif nomenclature
    :return: Motif, list(dict), int, int - motif, result lines, input length, length of filtered intput
    """
    # build motif class
    name = None if 'name' not in df.columns or df.iloc[0]['name'] in ['None', ''] else df.iloc[0]['name']
    motif_class = Motif(motif_str, name)

    # filter/cut those with low quality
    input_len = len(df)
    # TODO: extract to function?
    if args.skip_quality_under > 0 and 'quality' in df.columns:
        filtered_df = df[df.apply(
            lambda row: has_good_quality(row, args.skip_quality_under, 1, len(motif_class.get_modules()) - 2),
            axis=1
        )]
        log_str("Kept {:4d}/{:4d} ({:5.1f}%) reads for {}".format(
            len(filtered_df), len(df), len(filtered_df) / len(df) * 100.0, motif_class.dir_name()
        ))
        df = filtered_df
    # TODO: extract to function?
    if args.cut_quality_under > 0 and 'quality' in df.columns:
        cut_df: pd.DataFrame = df.apply(
            lambda row: cut_low_quality(row, args.cut_quality_under),
            result_type="expand",  # forcing DataFrame output, instead of Series
            axis=1
        )
        kept_bases = cut_df['read'].str.len().sum()
        all_bases = df['read'].str.len().sum()
        log_str("Cut {:4d}/{:4d} ({:5.1f}%) bases for {}".format(
            (all_bases - kept_bases), all_bases, ((all_bases - kept_bases) / all_bases * 100.0),
            motif_class.dir_name()
        ))
        df = cut_df

    filtered_len = len(df)
    if filtered_len == 0:
        # TODO: maybe it would be nice to warn user?
        return motif_class, [], input_len, filtered_len

    # https://github.com/pandas-dev/pandas-stubs/issues/1001
    # create annotations from rows
    annotations: list[Annotation] = []
    for _, row in df.iterrows():
        annotations.append(Annotation(
            row['read_id'], row['mate_order'], row['read'], row['reference'],
            row['modules'], row['log_likelihood'], motif_class
        ))

    # create annotation pairs from annotations
    annotation_pairs = annotations_to_pairs(annotations)

    # deduplicate?
    if args.deduplicate:
        annotation_pairs, duplicates = remove_pcr_duplicates(annotation_pairs)

    # infer read distribution
    read_distribution = np.bincount([len(ann.read_seq) for ann in annotations], minlength=100)

    # create report for each repeating module
    result_lines: list[dict] = []
    repeating_modules = motif_class.get_repeating_modules()
    postfilter_class = PostFilter(args)
    for i, (module_number, _, _) in enumerate(repeating_modules):
        prev_module = None if i == 0 else repeating_modules[i - 1]
        file_pcolor, file_output = None, None
        if args.verbose:
            motif_dir = f'{args.output_dir}/{motif_class.dir_name()}'
            os.makedirs(motif_dir, exist_ok=True)
            file_pcolor = f'{motif_dir}/pcolor_{module_number}'
            file_output = f'{motif_dir}/allcall_{module_number}.txt'

        annotations1: list[Annotation]
        annotations1 = pairs_to_annotations_pick(annotation_pairs, module_number)
        result = too_complex_function_inferring_and_creating_reports(
            module_number=module_number,
            postfilter_class=postfilter_class, motif_dir=motif_dir,
            motif_class=motif_class, read_distribution=read_distribution,
            motif_str=motif_str, result_lines=result_lines, prev_module=prev_module,
            args=args, file_pcolor=file_pcolor, file_output=file_output, annotations1=annotations1
        )
        postfilter_counts, postfilter_counts_phasing, \
            phasing, supp_reads, confidence, predicted, result_lines, last_num1 = result

        if args.verbose:
            report(
                args, motif_dir, motif_class, module_number, last_num1, prev_module,
                postfilter_counts[0], postfilter_counts[1], postfilter_counts[2],
                postfilter_counts_phasing[0], postfilter_counts_phasing[1], postfilter_counts_phasing[2],
                phasing, supp_reads, confidence, predicted
            )

    if args.verbose:
        # generate nomenclatures for all modules:
        for i, (seq, reps) in enumerate(motif_class.get_modules()):
            # write files if needed
            if reps == 1:
                # pick annotations from pairs if needed
                annotations = pairs_to_annotations_pick(annotation_pairs, i)

                # setup post filtering - no primers, insufficient quality, ...
                qual_annot, _ = postfilter_class.get_filtered(annotations, i, both_primers=True)

                # gather and write nomenclatures
                if len(qual_annot) > 0:
                    write_histogram_nomenclature(
                        f'{motif_dir}/nomenclatures_{i}.txt', qual_annot, index_rep=i
                    )
        # try to get the overall nomenclature:
        annotations = pairs_to_annotations_pick(annotation_pairs, None)
        for module_number, _, _ in repeating_modules:
            annotations, _ = postfilter_class.get_filtered(annotations, module_number, both_primers=True)
        write_histogram_nomenclature(f'{motif_dir}/nomenclature.txt', annotations)

    # return motif name in case it was processed normally
    return motif_class, result_lines, input_len, filtered_len


def process_group_tuple(
    x: tuple[argparse.Namespace, pd.DataFrame, str]
) -> tuple[Motif, list[dict], int, int]:
    """
    Wrapper for process_group() to use in parallelization (pool.imap).
    """
    return process_group(x[0], x[1], x[2])


def generate_groups_gzipped(
    input_stream: TextIO, column_name: str = 'motif', chunk_size: int = 1000000
) -> Iterator[pd.DataFrame]:
    """
    Generate sub-parts of the input table according to "column_name".
    Able to process even large files. Gzipped version
    :param input_stream: TextIO - input stream
    :param column_name: str - column name for grouping of the table
    :param chunk_size: int - chunk size for table processing
    :return: Iterator[pd.DataFrame] - sub-parts of the input table
    """
    with gzip.GzipFile(fileobj=input_stream.buffer, mode='r') as gz_file:
        # convert to text stream
        text_stream = io.TextIOWrapper(gz_file)

        for g in generate_groups(text_stream, column_name, chunk_size):
            yield g


def generate_groups(
    input_stream: TextIO, column_name: str = 'motif', chunk_size: int = 1000000
) -> Iterator[pd.DataFrame]:
    """
    Generate sub-parts of the input table according to "column_name".
    Able to process even large files.
    :param input_stream: TextIO - input stream
    :param column_name: str - column name for grouping of the table
    :param chunk_size: int - chunk size for table processing
    :return: Iterator[pd.DataFrame] - sub-parts of the input table
    """
    # initialize reading
    current_group_data = pd.DataFrame()

    # read the output of remaSTR into annotations
    for chunk in pd.read_csv(
        input_stream, sep='\t', chunksize=chunk_size, iterator=True, quoting=csv.QUOTE_NONE
    ):

        # identify the unique groups in the chunk
        unique_groups = chunk[column_name].unique()

        # got through the groups
        for group in unique_groups:

            # filter rows for the current group from the chunk
            group_data = chunk[chunk[column_name] == group]

            # if this group is a continuation of the previous group from the last chunk
            if not current_group_data.empty and current_group_data[column_name].iloc[0] == group:
                current_group_data = pd.concat([current_group_data, group_data])
                continue  # Move to the next group in the current chunk

            # if there's data in current_group_data, process and empty it
            if not current_group_data.empty:
                yield current_group_data
                current_group_data = pd.DataFrame()  # reset

            # check if the group is at end the chunk (probably continues into next one)
            if group == unique_groups[-1]:
                current_group_data = group_data
            else:
                # process the group as normally
                yield group_data

    # process the last group
    if not current_group_data.empty:
        yield current_group_data


def consume_iterator(
    results_iterator: Iterable[tuple[Motif, list[dict], int, int]]
) -> tuple[list[Motif], pd.DataFrame, int, int]:
    """
    Consume iterator of results.
    :param results_iterator: generator - motif and its corresponding results of modules
    :return: list[Motif], pd.DataFrame, int - motifs in list and table of all results, input length
    """
    # consume iterator of results
    all_motifs: list[Motif] = []
    all_result_lines: list[dict[Any, Any]] = []
    all_input_len = 0
    all_filtered_len = 0
    for i, (motif, rls, input_l, filtered_l) in enumerate(results_iterator):
        if filtered_l > 0:
            # append data
            all_motifs.append(motif)
            all_result_lines.extend(rls)
            all_input_len += input_l
            all_filtered_len += filtered_l

        # report progress
        # if args.progress > 0 and (i + 1) % args.progress == 0:
        #     report.log_str(f'Progress: {i + 1:10d} motifs done. ({datetime.now():%Y-%m-%d %H:%M:%S})')

    df = pd.DataFrame.from_records(all_result_lines)
    df.sort_values(by=['motif_name'], kind='stable')
    return (all_motifs, df, all_input_len, all_filtered_len)


def normalize_ref_alt(ref: str, alt: str) -> tuple[str, str]:
    suffix = 0
    it = zip(reversed(ref), reversed(alt))
    for _ in range(min(len(ref), len(alt)) - 1):
        (a, b) = next(it)
        if a == b:
            suffix += 1

    return (ref[:len(ref) - suffix], alt[:len(alt) - suffix])


def make_vcf_line(
    chrom: str, pos: str, unit: str, ref_copies: str, alt_copies: str, genotype: str,
    lines: list[str]
):
    if ref_copies == alt_copies:
        return  # they are the same, there is no variant
    if "|" in alt_copies:
        return  # skip over phased variants, TODO: later

    ref_seq = unit * int(ref_copies)
    if alt_copies == "B":
        info = f"REF={ref_copies};RU={unit};BG;END={pos}"
        lines.append("{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\n".format(
            chrom, pos, ".", ref_seq, "<BG>", ".", "PASS", info, "GT", genotype
        ))
    elif alt_copies == "E":
        info = f"REF={ref_copies};RU={unit};EXP;END={pos}"
        lines.append("{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\n".format(
            chrom, pos, ".", ref_seq[0], "<EXP>", ".", "PASS", info, "GT", genotype
        ))
    else:
        alt_seq = unit * int(alt_copies)

        svlen = len(alt_seq) - len(ref_seq)
        svtype = "INS" if svlen > 0 else "DEL"
        (ref_seq, alt_seq) = normalize_ref_alt(ref_seq, alt_seq)

        info = f"REF={ref_copies};RU={unit};SVLEN={svlen};SVTYPE={svtype}"
        lines.append("{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\n".format(
            chrom, pos, ".", ref_seq, alt_seq, ".", "PASS", info, "GT", genotype
        ))


def write_vcf(df: pd.DataFrame, out: str) -> None:
    lines = []
    lines.append('##fileformat=VCFv4.1\n')
    lines.append('##ALT=<ID=BG,Description="Background">\n')
    lines.append('##ALT=<ID=EXP,Description="Expansion of unknown (large) size">\n')
    lines.append('##FILTER=<ID=PASS,Description="All filters passed">\n')
    lines.append('##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">\n')
    lines.append('##INFO=<ID=END,Number=1,Type=Integer,Description="End position of the variant">\n')
    lines.append('##INFO=<ID=BG,Number=0,Type=Flag,Description="Background variant">\n')
    lines.append('##INFO=<ID=EXP,Number=0,Type=Flag,Description="Expansion variant">\n')
    lines.append('##INFO=<ID=REF,Number=1,Type=Integer,Description="Reference copy number">\n')
    lines.append('##INFO=<ID=RU,Number=1,Type=String,Description="Repeat unit in ref orientation">\n')
    lines.append('##INFO=<ID=SVLEN,Number=1,Type=Integer,Description="Alt length - Ref length">\n')
    lines.append('##INFO=<ID=SVTYPE,Number=1,Type=String,Description="Type of structural variant">\n')
    lines.append('#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tsample\n')

    records: list[str] = []
    for i, row in df.iterrows():
        m1 = re.match(r"([A-Z]+)\[([0-9]+)\]", row["motif_sequence"])
        if m1 is None:
            print(f"{row['motif_sequence']} returned None")
            continue
        unit, copies = m1.groups()
        allele1 = str(row["allele1"])
        allele2 = str(row["allele2"])

        if allele1 == allele2:
            make_vcf_line(row["chromosome"], row["start"], unit, copies, str(row["allele1"]), "1/1", records)
        else:
            make_vcf_line(row["chromosome"], row["start"], unit, copies, str(row["allele1"]), "1/.", records)
            make_vcf_line(row["chromosome"], row["start"], unit, copies, str(row["allele2"]), "./1", records)

    records.sort(key=lambda x: chr_and_pos(x))
    with open(f"{out}/variants.vcf", "w") as f:
        f.writelines(lines + records)


def chr_and_pos(line: str) -> tuple[int, int]:
    m1 = re.match(r"(chr[0-9XY]+)\t([0-9]+)\t.*", line)
    if m1 is None:
        raise ValueError(f"got {line}")
    chrom, pos = m1.groups()
    chrom2: int = {
        "chr1":  1,  "chr2":  2,  "chr3":  3,  "chr4":  4,  "chr5": 5,   "chr6": 6,
        "chr7":  7,  "chr8":  8,  "chr9":  9,  "chr10": 10, "chr11": 11, "chr12": 12,
        "chr13": 13, "chr14": 14, "chr15": 15, "chr16": 16, "chr17": 17, "chr18": 18,
        "chr19": 19, "chr20": 20, "chr21": 21, "chr22": 22, "chrX": 23,  "chrY": 24
    }[chrom]
    return (chrom2, int(pos))


def main() -> None:
    start_time = datetime.now()
    args = load_arguments()

    # initialize logging module
    configure_logger(f'{args.output_dir}/dante.log')
    log_str('DANTE_remaSTR = "Da Amazing NucleoTide Exposer" (remastered)')
    log_str(f'DANTE_remaSTR Starting : {start_time:%Y-%m-%d %H:%M:%S}')

    # process the input
    motif_column_name = 'motif'
    if args.input_gzipped:
        groups_iterator = generate_groups_gzipped(sys.stdin, motif_column_name)
    else:
        groups_iterator = generate_groups(sys.stdin, motif_column_name)

    stop_motif = args.start_motif + args.max_motifs if args.max_motifs is not None else None
    groups_iterator = itertools.islice(groups_iterator, args.start_motif, stop_motif)

    # create iterator of results
    all_inputs = ((args, motif_table, motif_table[motif_column_name].iloc[0])
                  for motif_table in groups_iterator)

    iterator: Iterable
    if args.processes > 1:
        with multiprocessing.Pool(args.processes) as pool:
            iterator = pool.imap(process_group_tuple, all_inputs, chunksize=5)
            output = consume_iterator(iterator)
    else:
        iterator = (process_group_tuple(inputs) for inputs in all_inputs)
        output = consume_iterator(iterator)

    all_motifs, rl_df, input_len, filtered_len = output

    # summary of the filtered reads
    log_str(f'Kept {filtered_len:4d}/{input_len:4d} ({filtered_len / input_len * 100.0:5.1f}%) reads.')

    #  write the dataframe to stdout
    out_file = args.output_dir + "/variants.tsv"
    log_str(f'Writing output to {out_file}: {datetime.now():%Y-%m-%d %H:%M:%S}')
    rl_df.to_csv(out_file, sep='\t')

    write_vcf(rl_df, args.output_dir)

    # generate report and output files for the whole run
    # generate report.html and alignments.html, but not fastas
    if args.verbose:
        post_filter = PostFilter(args)
        write_report(
            sorted(all_motifs), rl_df, post_filter,
            args.output_dir, args.nomenclatures
        )

    # print the time of the end
    end_time = datetime.now()
    log_str(f'DANTE_remaSTR Stopping : {end_time:%Y-%m-%d %H:%M:%S}')
    log_str(f'Total time of run      : {end_time - start_time}')


if __name__ == '__main__':
    # this is good practice so the variables do not pollute the global scope
    main()
