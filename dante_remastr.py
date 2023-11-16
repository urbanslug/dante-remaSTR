import argparse
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

MOTIF_PRINT_LEN = 40
result_line_templ = ('{motif_name}\t{motif_nom}\t{motif_seq}\t{motif_chrom}\t{motif_start}\t{motif_end}\t{allele1}\t{allele2}\t{confidence}\t'
                     '{confidence1}\t{confidence2}\t{qual_num}\t{primer_num}\t{filt_num}\t{conf_back}\t{conf_back_all}\t{conf_ext}\t{conf_ext_all}')


def shorten_str(string: str, max_length: int = MOTIF_PRINT_LEN, ellipsis_str: str = '...') -> str:
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


def generate_result_line(motif_class: Motif, predicted: tuple[str, str], confidence: tuple[float | str, ...], qual_num: int,
                         primer_num: int, filt_num: int, module_number: int, second_module_number: int = None, ):
    """
    Generate result line from the template string.
    :param motif_class: Motif - motif class
    :param predicted: tuple[str, str] - predicted alleles (number or 'B'/'E')
    :param confidence: tuple[7x float/str] - confidences of prediction
    :param qual_num: int - number of reads with both primers
    :param primer_num: int - number of reads with exactly one primer
    :param filt_num: int - number of filtered out reads (no primers, many errors, ...)
    :param module_number: int - module number in motif
    :param second_module_number: int/None - second module number in motif
    :return: str - tsv line with result
    """
    start, end = motif_class.get_location_subpart(module_number)
    motif_seq = motif_class.module_str(module_number)
    if second_module_number is not None:
        _, end = motif_class.get_location_subpart(second_module_number)
        motif_seq = ','.join([motif_class.module_str(i) for i in range(module_number, second_module_number + 1)])
    return result_line_templ.format(motif_name=motif_class.name, motif_nom=motif_class.motif, motif_seq=motif_seq,
                                    motif_chrom=motif_class.chrom, motif_start=start, motif_end=end, allele1=predicted[0], allele2=predicted[1],
                                    confidence=confidence[0], confidence1=confidence[1], confidence2=confidence[2], qual_num=qual_num,
                                    primer_num=primer_num, filt_num=filt_num, conf_back=confidence[3] if len(confidence) > 3 else '---',
                                    conf_back_all=confidence[4] if len(confidence) > 4 else '---',
                                    conf_ext=confidence[5] if len(confidence) > 5 else '---',
                                    conf_ext_all=confidence[6] if len(confidence) > 6 else '---')


def generate_result_header():
    """
    Generate result header from the template string.
    :return: str - tsv header
    """
    return result_line_templ.format(motif_name='motif_name', motif_nom='motif_nomenclature', motif_seq='motif_sequence', motif_chrom='chromosome',
                                    motif_start='start', motif_end='end', allele1='allele1', allele2='allele2', confidence='confidence',
                                    confidence1='conf_allele1', confidence2='conf_allele2', qual_num='quality_reads', primer_num='one_primer_reads',
                                    filt_num='filtered_reads', conf_back='conf_background', conf_back_all='conf_background_all',
                                    conf_ext='conf_extended', conf_ext_all='conf_extended_all')


def process_group(args: argparse.Namespace, df: pd.DataFrame, motif_str: str) -> tuple[Motif, list[str]]:
    """
    Process the group as pandas Dataframe. Return motif name if processed correctly or None otherwise.
    :param args: argparse.Namespace - namespace of program arguments
    :param df: pd.DataFrame - pandas DataFrame with information about annotated reads for a single motif to process
    :param motif_str: str - motif nomenclature
    :return: Motif, list(str) - motif, result lines
    """
    # build motif class
    name = None if 'name' not in df.columns or df.iloc[0]['name'] in ['None', ''] else df.iloc[0]['name']
    motif_class = Motif(motif_str, name)

    # setup motif_directory
    motif_dir = f'{args.output_dir}/{motif_class.dir_name()}'

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
            report.write_all(qual_annot, primer_annot, filt_annot, motif_dir, motif_class, module_number, zip_it=args.gzip_outputs)

        # run inference - this takes most of the time (for no --verbose)
        file_pcolor = f'{motif_dir}/pcolor_{module_number}' if args.verbose else None
        file_output = f'{motif_dir}/allcall_{module_number}.txt' if args.verbose else None
        inference_class = inference.Inference(read_distribution, args.param_file, str_rep=args.min_rep_cnt, minl_primer1=args.min_flank_len,
                                              minl_primer2=args.min_flank_len, minl_str=args.min_rep_len)
        predicted, confidence = inference_class.genotype(qual_annot, primer_annot, module_number, file_pcolor, file_output, motif_str)

        # write the alignments
        if confidence is not None and args.verbose:
            conf, c1, c2, _, _, _, _ = confidence
            a1 = int(predicted[0]) if isinstance(predicted[0], int) else None
            a2 = int(predicted[1]) if isinstance(predicted[1], int) else None
            if a1 is not None and a1 > 0:
                report.write_alignment(f'{motif_dir}/alignment_{module_number}_a{a1}.fasta', qual_annot, module_number, a1, zip_it=args.gzip_outputs)
            if a2 is not None and a2 != a1 and a2 != 0:
                report.write_alignment(f'{motif_dir}/alignment_{module_number}_a{a2}.fasta', qual_annot, module_number, a2, zip_it=args.gzip_outputs)

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
                                 zip_it=args.gzip_outputs)

            # infer phasing
            phasing, supp_reads = inference.phase(both_good_annot, last_num, module_number)

            # write phasing into a file
            if args.verbose:
                inference.save_phasing(f'{motif_dir}/phasing_{last_num}_{module_number}.txt', phasing, supp_reads)

            # append to the result line TODO remove from here? and add to html exports.
            result_lines.append(
                generate_result_line(motif_class, phasing, supp_reads, len(both_good_annot), len(one_good_annot), len(none_good_annot), last_num,
                                     module_number))

        # append to the result line
        result_lines.append(
            generate_result_line(motif_class, predicted, confidence, len(qual_annot), len(primer_annot), len(filt_annot), module_number))

    # try to get the overall nomenclature:
    if args.verbose:
        annotations = annotation.pairs_to_annotations_pick(annotation_pairs, None)
        for module_number, _, _ in motif_class.get_repeating_modules():
            postfilter_class = PostFilter(args)
            annotations, _ = postfilter_class.get_filtered(annotations, module_number, both_primers=True)
        report.write_histogram_nomenclature(f'{motif_dir}/nomenclature.txt', annotations)

    # return motif name in case it was processed normally
    return motif_class, result_lines


def process_group_tuple(x: tuple[argparse.Namespace, pd.DataFrame, str]) -> tuple[Motif, list[str]]:
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
    for chunk in pd.read_csv(input_stream, sep='\t', chunksize=chunk_size, iterator=True):

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


def report_group(result_lines: list[str], progress: int, num_processed: int) -> None:
    """
    Report results of a group.
    :param result_lines: list(str) - result lines to report
    :param progress: int - how many lines to report progress
    :param num_processed: int - how many were processed already
    """
    for result_line in result_lines:
        report.log_str(result_line, stdout_too=sys.stdout)
    if progress > 0 and num_processed % progress == 0:
        report.log_str(f'Progress: {num_processed:10d} motifs done. ({datetime.now():%Y-%m-%d %H:%M:%S})')


if __name__ == '__main__':
    # save the time of the start
    start_time = datetime.now()

    # load arguments
    args = arguments.load_arguments()

    # initialize logging module
    report.configure_logger(f'{args.output_dir}/dante.log')
    report.log_str('DANTE_remaSTR = "Da Amazing NucleoTide Exposer" (remastered)')
    report.log_str(f'DANTE_remaSTR Starting : {start_time:%Y-%m-%d %H:%M:%S}')

    # log header, data will be logged in the process_group
    report.log_str(generate_result_header(), stdout_too=sys.stdout)

    # process the input
    motif_column_name = 'motif'
    processed_motifs = []
    groups_iterator = generate_groups_gzipped(sys.stdin, motif_column_name) if args.input_gzipped else generate_groups(sys.stdin, motif_column_name)
    groups_iterator = itertools.islice(groups_iterator, args.start_motif, args.start_motif + args.max_motifs if args.max_motifs is not None else None)
    all_inputs = ((args, motif_table, motif_table[motif_column_name].iloc[0]) for motif_table in groups_iterator)
    if args.processes > 1:
        with multiprocessing.Pool(args.processes) as pool:
            for motif, rls in pool.imap(process_group_tuple, all_inputs, chunksize=100):
                processed_motifs.append(motif)
                report_group(rls, args.progress, len(processed_motifs))
    else:
        for motif, rls in (process_group(*inputs) for inputs in all_inputs):
            processed_motifs.append(motif)
            report_group(rls, args.progress, len(processed_motifs))

    # generate report and output files for the whole run
    if args.verbose:
        post_filter = PostFilter(args)
        report.write_report(processed_motifs, post_filter, args.output_dir, args.nomenclatures)

    # print the time of the end
    end_time = datetime.now()
    report.log_str(f'DANTE_remaSTR Stopping : {end_time:%Y-%m-%d %H:%M:%S}')
    report.log_str(f'Total time of run      : {end_time - start_time}')
