import argparse
import csv
import gzip
import io
import sys
import typing
import itertools
import multiprocessing
from datetime import datetime

import numpy as np
import pandas as pd

import src.inference as inference
import src.arguments as arguments
import src.report as report
import src.annotation as annotation
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


def errors_per_read(errors: list[tuple[int, int, int]], relative: bool = False) -> tuple[float | str, float | str]:
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
        return (np.mean([indels / float(length) for indels, _, length in errors]) * mean_length,
                np.mean([mismatches / float(length) for _, mismatches, length in errors]) * mean_length)
    else:
        return np.mean([indels for indels, _, _ in errors]), np.mean([mismatches for _, mismatches, _ in errors])


def generate_result_line(motif_class: Motif, predicted: tuple[str, str], confidence: tuple[float | str, ...], qual_num: int, primer_num: int,
                         filt_num: int, module_number: int, qual_annot: list[annotation.Annotation] = None, second_module_number: int = None) -> dict:
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
        motif_seq = ','.join([motif_class.module_str(i) for i in range(module_number, second_module_number + 1)])

    # get info about errors and number of reads from quality annotations if provided
    reads_a1 = reads_a2 = '---'
    indels_rel = mismatches_rel = '---'
    indels_rel1 = mismatches_rel1 = '---'
    indels_rel2 = mismatches_rel2 = '---'
    if qual_annot is not None:
        # get info about number of reads
        a1 = int(predicted[0]) if isinstance(predicted[0], int) else None
        a2 = int(predicted[1]) if isinstance(predicted[1], int) else None
        reads_a1 = 0 if a1 is None else len([a for a in qual_annot if a.module_repetitions[module_number] == a1])
        reads_a2 = 0 if a2 is None else len([a for a in qual_annot if a.module_repetitions[module_number] == a2])

        # get info about errors
        errors = [a.get_module_errors(module_number) for a in qual_annot]
        errors_a1 = [a.get_module_errors(module_number) for a in qual_annot if a.module_repetitions[module_number] == a1]
        errors_a2 = [a.get_module_errors(module_number) for a in qual_annot if a.module_repetitions[module_number] == a2]
        assert len([l for i, m, l in errors if l == 0]) == 0

        # extract error metrics
        indels_rel, mismatches_rel = errors_per_read(errors, relative=True)
        indels_rel1, mismatches_rel1 = errors_per_read(errors_a1, relative=True)
        indels_rel2, mismatches_rel2 = errors_per_read(errors_a2, relative=True)

    # return dictionary
    return {'motif_name': motif_class.name, 'motif_nomenclature': motif_class.motif, 'motif_sequence': motif_seq, 'chromosome': motif_class.chrom,
            'start': start, 'end': end, 'allele1': predicted[0], 'allele2': predicted[1], 'confidence': confidence[0], 'conf_allele1': confidence[1],
            'conf_allele2': confidence[2], 'reads_a1': reads_a1, 'reads_a2': reads_a2, 'indels': indels_rel, 'mismatches': mismatches_rel,
            'indels_a1': indels_rel1, 'indels_a2': indels_rel2, 'mismatches_a1': mismatches_rel1, 'mismatches_a2': mismatches_rel2,
            'quality_reads': qual_num, 'one_primer_reads': primer_num, 'filtered_reads': filt_num,
            'conf_background': confidence[3] if len(confidence) > 3 else '---',
            'conf_background_all': confidence[4] if len(confidence) > 4 else '---',
            'conf_extended': confidence[5] if len(confidence) > 5 else '---', 'conf_extended_all': confidence[6] if len(confidence) > 6 else '---',
            'repetition_index': module_number if second_module_number is None else f'{module_number}_{second_module_number}'}


def process_group(args: argparse.Namespace, df: pd.DataFrame, motif_str: str) -> tuple[Motif, list[dict], int, int]:
    """
    Process the group as pandas Dataframe. Return motif name if processed correctly or None otherwise.
    :param args: argparse.Namespace - namespace of program arguments
    :param df: pd.DataFrame - pandas DataFrame with information about annotated reads for a single motif to process
    :param motif_str: str - motif nomenclature
    :return: Motif, list(dict), int - motif, result lines, input length, length of filtered intput
    """
    # build motif class
    name = None if 'name' not in df.columns or df.iloc[0]['name'] in ['None', ''] else df.iloc[0]['name']
    motif_class = Motif(motif_str, name)

    # setup motif_directory
    motif_dir = f'{args.output_dir}/{motif_class.dir_name()}'

    # filter/cut those with low quality
    input_len = len(df)
    if args.skip_quality_under > 0 and 'quality' in df.columns:
        filtered_df = df[df.apply(lambda row: has_good_quality(row, args.skip_quality_under, 1, len(motif_class.get_modules()) - 2), axis=1)]
        report.log_str(f'Kept {len(filtered_df):4d}/{len(df):4d} ({len(filtered_df) / len(df) * 100.0:5.1f}%) reads for {motif_class.dir_name()}')
        df = filtered_df
    if args.cut_quality_under > 0 and 'quality' in df.columns:
        cut_df = df.apply(lambda row: cut_low_quality(row, args.cut_quality_under), axis=1)
        kept_bases = cut_df['read'].str.len().sum()
        all_bases = df['read'].str.len().sum()
        report.log_str(f'Cut {all_bases - kept_bases:4d}/{all_bases:4d} ({(all_bases - kept_bases) / all_bases * 100.0:5.1f}%) bases for '
                       f'{motif_class.dir_name()}')
        df = cut_df
    filtered_len = len(df)

    # create annotations from rows
    annotations = df.apply(
        lambda row: annotation.Annotation(row['read_id'], row['mate_order'], row['read'], row['reference'], row['modules'],
                                          row['log_likelihood'], motif_class), axis=1)

    # create annotation pairs from annotations
    annotation_pairs = annotation.annotations_to_pairs(annotations)

    # deduplicate?
    if args.deduplicate:
        annotation_pairs, duplicates = annotation.remove_pcr_duplicates(annotation_pairs)

    # infer read distribution
    read_distribution = np.bincount([len(ann.read_seq) for ann in annotations], minlength=100)

    # create report for each repeating module
    result_lines = []
    repeating_modules = motif_class.get_repeating_modules()
    for i, (module_number, _, _) in enumerate(repeating_modules):

        # pick annotations from pairs if needed
        annotations = annotation.pairs_to_annotations_pick(annotation_pairs, module_number)

        # setup post filtering - no primers, insufficient quality, ...
        postfilter_class = PostFilter(args)
        qual_annot, filt_annot = postfilter_class.get_filtered(annotations, module_number, both_primers=True)
        primer_annot, filt_annot = postfilter_class.get_filtered(filt_annot, module_number, both_primers=False)

        # write files if needed
        if args.verbose:
            report.write_all(qual_annot, primer_annot, filt_annot, motif_dir, motif_class, module_number,
                             zip_it=args.gzip_outputs, cutoff_alignments=args.cutoff_alignments)

        # run inference - this takes most of the time (for no --verbose)
        file_pcolor = f'{motif_dir}/pcolor_{module_number}' if args.verbose else None
        file_output = f'{motif_dir}/allcall_{module_number}.txt' if args.verbose else None
        inference_class = inference.Inference(read_distribution, args.param_file, str_rep=args.min_rep_cnt, minl_primer1=args.min_flank_len,
                                              minl_primer2=args.min_flank_len, minl_str=args.min_rep_len)
        monoallelic = args.male and report.chrom_from_string(motif_class.chrom) in [report.ChromEnum.X, report.ChromEnum.Y]
        predicted, confidence = inference_class.genotype(qual_annot, primer_annot, module_number, file_pcolor, file_output, motif_str, monoallelic)

        # get number of precise alignments for each allele
        a1 = int(predicted[0]) if isinstance(predicted[0], int) else None
        a2 = int(predicted[1]) if isinstance(predicted[1], int) else None

        # write the alignments
        if confidence is not None and args.verbose:
            if a1 is not None and a1 > 0:
                report.write_alignment(f'{motif_dir}/alignment_{module_number}_a{a1}.fasta', qual_annot, module_number, a1,
                                       zip_it=args.gzip_outputs, cutoff_after=args.cutoff_alignments)
            if a2 is not None and a2 != a1 and a2 != 0:
                report.write_alignment(f'{motif_dir}/alignment_{module_number}_a{a2}.fasta', qual_annot, module_number, a2,
                                       zip_it=args.gzip_outputs, cutoff_after=args.cutoff_alignments)

        # infer phasing (if we are not on the first repeating module)
        if i != 0:
            # get the last module number
            last_num = repeating_modules[i - 1][0]

            # post filtering
            both_good_annot, filtered_annot = postfilter_class.get_filtered_list(annotations, [last_num, module_number], both_primers=[True, True])
            left_good_annot, left_bad_annot = postfilter_class.get_filtered_list(filtered_annot, [last_num, module_number],
                                                                                 both_primers=[False, True])
            right_good_annot, none_good_annot = postfilter_class.get_filtered_list(left_bad_annot, [last_num, module_number],
                                                                                   both_primers=[True, False])
            one_good_annot = left_good_annot + right_good_annot

            # write files
            if args.verbose:
                report.write_all(both_good_annot, one_good_annot, none_good_annot, motif_dir, motif_class, last_num, module_number,
                                 zip_it=args.gzip_outputs, cutoff_alignments=args.cutoff_alignments)

            # infer phasing
            phasing, supp_reads = inference.phase(both_good_annot, last_num, module_number)

            # write phasing into a file
            if args.verbose:
                inference.save_phasing(f'{motif_dir}/phasing_{last_num}_{module_number}.txt', phasing, supp_reads)

            # append to the result line
            result_lines.append(
                generate_result_line(motif_class, phasing, supp_reads, len(both_good_annot), len(one_good_annot), len(none_good_annot), last_num,
                                     second_module_number=module_number))

        # append to the result line
        result_lines.append(
            generate_result_line(motif_class, predicted, confidence, len(qual_annot), len(primer_annot), len(filt_annot), module_number,
                                 qual_annot=qual_annot))

    # try to get the overall nomenclature:
    if args.verbose:
        annotations = annotation.pairs_to_annotations_pick(annotation_pairs, None)
        for module_number, _, _ in repeating_modules:
            postfilter_class = PostFilter(args)
            annotations, _ = postfilter_class.get_filtered(annotations, module_number, both_primers=True)
        report.write_histogram_nomenclature(f'{motif_dir}/nomenclature.txt', annotations)

    # return motif name in case it was processed normally
    return motif_class, result_lines, input_len, filtered_len


def process_group_tuple(x: tuple[argparse.Namespace, pd.DataFrame, str]) -> tuple[Motif, list[dict], int, int]:
    """
    Wrapper for process_group() to use in parallelization (pool.imap).
    """
    return process_group(x[0], x[1], x[2])


def generate_groups_gzipped(input_stream: typing.TextIO, column_name: str = 'motif', chunk_size: int = 1000000) -> typing.Iterator[pd.DataFrame]:
    """
    Generate sub-parts of the input table according to "column_name". Able to process even large files. Gzipped version
    :param input_stream: TextIO - input stream
    :param column_name: str - column name for grouping of the table
    :param chunk_size: int - chunk size for table processing
    :return: Iterator[pd.DataFrame] - sub-parts of the input table
    """
    with gzip.GzipFile(fileobj=input_stream.buffer, mode='r') as gz_file:  # type: typing.IO[bytes]
        # convert to text stream
        text_stream = io.TextIOWrapper(gz_file)

        for g in generate_groups(text_stream, column_name, chunk_size):
            yield g


def generate_groups(input_stream: typing.TextIO, column_name: str = 'motif', chunk_size: int = 1000000) -> typing.Iterator[pd.DataFrame]:
    """
    Generate sub-parts of the input table according to "column_name". Able to process even large files.
    :param input_stream: TextIO - input stream
    :param column_name: str - column name for grouping of the table
    :param chunk_size: int - chunk size for table processing
    :return: Iterator[pd.DataFrame] - sub-parts of the input table
    """
    # initialize reading
    current_group_data = pd.DataFrame()

    # read the output of remaSTR into annotations
    for chunk in pd.read_csv(input_stream, sep='\t', chunksize=chunk_size, iterator=True, quoting=csv.QUOTE_NONE):

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


def consume_iterator(results_iterator: typing.Generator[tuple[Motif, list[dict], int, int], any, None]) -> tuple[list[Motif], pd.DataFrame, int, int]:
    """
    Consume iterator of results.
    :param results_iterator: generator - motif and its corresponding results of modules
    :return: list[Motif], pd.DataFrame, int - motifs in list and table of all results, input length
    """
    # consume iterator of results
    all_result_lines = []
    all_motifs = []
    all_input_len = 0
    all_filtered_len = 0
    for i, (motif, rls, input_l, filtered_l) in enumerate(results_iterator):
        # append data
        all_motifs.append(motif)
        all_result_lines.extend(rls)
        all_input_len += input_l
        all_filtered_len += filtered_l

        # report progress
        if args.progress > 0 and (i + 1) % args.progress == 0:
            report.log_str(f'Progress: {i + 1:10d} motifs done. ({datetime.now():%Y-%m-%d %H:%M:%S})')

    return all_motifs, pd.DataFrame.from_dict(all_result_lines).sort_values(by=['motif_name'], kind='stable'), all_input_len, all_filtered_len


if __name__ == '__main__':
    # save the time of the start
    start_time = datetime.now()

    # load arguments
    args = arguments.load_arguments()

    # initialize logging module
    report.configure_logger(f'{args.output_dir}/dante.log')
    report.log_str('DANTE_remaSTR = "Da Amazing NucleoTide Exposer" (remastered)')
    report.log_str(f'DANTE_remaSTR Starting : {start_time:%Y-%m-%d %H:%M:%S}')

    # process the input
    motif_column_name = 'motif'
    groups_iterator = generate_groups_gzipped(sys.stdin, motif_column_name) if args.input_gzipped else generate_groups(sys.stdin, motif_column_name)
    groups_iterator = itertools.islice(groups_iterator, args.start_motif, args.start_motif + args.max_motifs if args.max_motifs is not None else None)
    all_inputs = ((args, motif_table, motif_table[motif_column_name].iloc[0]) for motif_table in groups_iterator)

    # create iterator of results
    if args.processes > 1:
        with multiprocessing.Pool(args.processes) as pool:
            all_motifs, rl_df, input_len, filtered_len = consume_iterator(pool.imap(process_group_tuple, all_inputs, chunksize=5))
    else:
        all_motifs, rl_df, input_len, filtered_len = consume_iterator((process_group(*inputs) for inputs in all_inputs))

    # summary of the filtered reads
    report.log_str(f'Kept {filtered_len:4d}/{input_len:4d} ({filtered_len / input_len * 100.0:5.1f}%) reads.')

    #  write the dataframe to stdout
    report.log_str(f'Writing output (stdout): {datetime.now():%Y-%m-%d %H:%M:%S}')
    rl_df.to_csv(sys.stdout, sep='\t')

    # generate report and output files for the whole run
    if args.verbose:
        post_filter = PostFilter(args)
        report.write_report(sorted(all_motifs), rl_df, post_filter, args.output_dir, args.nomenclatures)

    # print the time of the end
    end_time = datetime.now()
    report.log_str(f'DANTE_remaSTR Stopping : {end_time:%Y-%m-%d %H:%M:%S}')
    report.log_str(f'Total time of run      : {end_time - start_time}')
