import argparse
import gzip
import io
import sys
import typing
from datetime import datetime

import numpy as np
import pandas as pd

import all_call
import arguments
import postfilter
import report
import annotation
from annotation.Motif import Motif
from postfilter import PostFilter

MOTIF_PRINT_LEN = 40
result_line_template = ('{motif_name}\t{motif_seq}\t{motif_chrom}\t{motif_start}\t{motif_end}\t{allele1}\t{allele2}\t{confidence}\t{confidence1}\t'
                        '{confidence2}\t{qual_num}\t{primer_num}\t{filt_num}\t{conf_back}\t{conf_back_all}\t{conf_ext}\t{conf_ext_all}')


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


def generate_result_line(motif: Motif, module_number: int, predicted: tuple[str, str],
                         confidence: tuple[float, float, float, float, float, float, float], qual_num: int, primer_num: int, filt_num: int):
    seq, rep = motif[module_number]
    start, end = motif.get_location_subpart(module_number)
    return result_line_template.format(motif_name=motif.name, motif_seq=f'{seq}[{rep}]', motif_chrom=motif.chrom, motif_start=start,
                                       motif_end=end, allele1=predicted[0], allele2=predicted[1], confidence=confidence[0],
                                       confidence1=confidence[1], confidence2=confidence[2], qual_num=qual_num, primer_num=primer_num,
                                       filt_num=filt_num, conf_back=confidence[3], conf_back_all=confidence[4], conf_ext=confidence[5],
                                       conf_ext_all=confidence[6])


def generate_result_header():
    return result_line_template.format(motif_name='motif_name', motif_seq='motif_sequence', motif_chrom='chromosome', motif_start='start',
                                       motif_end='end', allele1='allele1', allele2='allele2', confidence='confidence', confidence1='conf_allele1',
                                       confidence2='conf_allele2', qual_num='quality_reads', primer_num='one_primer_reads', filt_num='filtered_reads',
                                       conf_back='conf_background', conf_back_all='conf_background_all', conf_ext='conf_extended',
                                       conf_ext_all='conf_extended_all')


def process_group(args: argparse.Namespace, df: pd.DataFrame, motif: str) -> typing.Optional[Motif]:
    """
    Process the group as pandas Dataframe. Return motif name if processed correctly or None otherwise.
    :param args: argparse.Namespace - namespace of program arguments
    :param df: pd.DataFrame - pandas DataFrame with information about annotated reads for a single motif to process
    :param motif: str - motif nomenclature
    :return: Motif/None - motif if processed correctly, None otherwise
    """
    # log header if we are at the beginning
    if process_group.launched == 0:
        report.log_str(generate_result_header(), stdout_too=sys.stdout)

    # count how many groups we have processed and end if enough/not already on start
    process_group.launched += 1
    if args.start_motif is not None and args.start_motif >= process_group.launched:
        return None
    if args.max_motifs is not None and args.max_motifs >= process_group.processed:
        return None
    process_group.processed += 1

    # build motif class
    motif_class = Motif(motif)

    # setup motif_directory
    motif_dir = f'{args.output_dir}/{motif_class.dir_name()}'
    # os.makedirs(motif_dir, exist_ok=True)

    # create annotations from rows
    annotations = df.apply(lambda row: annotation.Annotation(row['read_id'], row['read'], row['reference'], row['modules'], row['log_likelihood'],
                                                             motif_class), axis=1)

    # deduplicate?
    dedup_annot_pairs = []
    if args.deduplicate:
        # create annotation pairs from annotations
        annotation_pairs = annotation.annotations_to_pairs(annotations)

        # deduplicate
        dedup_annot_pairs, duplicates = annotation.remove_pcr_duplicates(annotation_pairs)

    # infer read distribution
    read_distribution = np.bincount([len(ann.read_seq) for ann in annotations], minlength=100)

    # create report for each repeating module
    for module_number, seq, _ in motif_class.get_repeating_modules():

        if args.deduplicate:
            # deduplicated annotations (for each pair we keep only one)
            annotations = annotation.pairs_to_annotations_pick(dedup_annot_pairs, module_number)

        # setup post filtering - no primers, insufficient quality, ...
        postfilter_class = postfilter.PostFilter(args)
        qual_annot, filt_annot = postfilter_class.get_filtered(annotations, module_number, both_primers=True)
        primer_annot, filt_annot = postfilter_class.get_filtered(filt_annot, module_number, both_primers=False)

        # write files if needed
        if args.verbose:
            report.write_all(qual_annot, primer_annot, filt_annot, motif_dir, module_number, zip_it=args.gzip_outputs)

        # run inference
        inference = all_call.Inference(read_distribution, args.param_file, str_rep=args.min_rep_cnt, minl_primer1=args.min_flank_len,
                                       minl_primer2=args.min_flank_len, minl_str=args.min_rep_len)
        file_pcolor = f'{motif_dir}/pcolor_{module_number}' if args.verbose else None
        file_output = f'{motif_dir}/allcall_{module_number}.txt' if args.verbose else None
        predicted, confidence = inference.all_call(qual_annot, primer_annot, module_number, file_pcolor, file_output, motif)

        # write the alignments
        if confidence is not None and args.verbose:
           conf, c1, c2, _, _, _, _ = confidence
           a1 = int(predicted[0]) if predicted[0].isdigit() else None
           a2 = int(predicted[1]) if predicted[1].isdigit() else None
           if a1 is not None and a1 > 0:
               report.write_alignment(f'{motif_dir}/alignment_{module_number}_a{a1}.fasta', qual_annot, module_number, allele=a1)
           if a2 is not None and a2 != a1 and a2 != 0:
               report.write_alignment(f'{motif_dir}/alignment_{module_number}_a{a2}.fasta', qual_annot, module_number, allele=a2)

        # print the output
        result_line = generate_result_line(motif_class, module_number, predicted, confidence, len(qual_annot), len(primer_annot), len(filt_annot))
        report.log_str(result_line, stdout_too=sys.stdout)

    # try to get the overall nomenclature:
    if args.verbose:
        if args.deduplicate:
            annotations = annotation.pairs_to_annotations_pick(dedup_annot_pairs, None)
        for module_number, seq, _ in motif_class.get_repeating_modules():
           postfilter_class = postfilter.PostFilter(args)
           qual_annot_all, _ = postfilter_class.get_filtered(annotations, module_number, both_primers=True)
        report.write_histogram_nomenclature(f'{motif_dir}/nomenclature.txt', annotations)

    # return motif name in case it was processed normally
    return motif_class

process_group.launched = 0
process_group.processed = 0


def generate_groups_gzipped(input_stream: typing.TextIO, column_name: str = 'motif', chunk_size: int = 1000000) -> typing.Iterator[pd.DataFrame]:
    """
    Generate sub-parts of the input table according to "column_name". Able to process even large files. Gzipped version
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
    processed_motifs = []
    for motif_table in generate_groups_gzipped(sys.stdin, motif_column_name) if args.input_gzipped else generate_groups(sys.stdin, motif_column_name):
        assert len(motif_table) > 0
        motif_name = motif_table[motif_column_name].iloc[0]
        processed_motif = process_group(args, motif_table, motif_name)
        if processed_motif is not None:
            processed_motifs.append(processed_motif)

    # generate report and output files for the whole run
    if args.verbose:
        post_filter = PostFilter(args)
        report.write_report(processed_motifs, post_filter, args.output_dir, args.nomenclatures)

    # print the time of the end
    end_time = datetime.now()
    report.log_str(f'DANTE_remaSTR Stopping : {end_time:%Y-%m-%d %H:%M:%S}')
    report.log_str(f'Total time of run      : {end_time - start_time}')
