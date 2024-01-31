import argparse
import glob
import os
import re
import textwrap

import pandas as pd


def load_arguments() -> argparse.Namespace:
    """
    Loads and parses arguments.
    :return: args - parsed arguments
    """
    parser = argparse.ArgumentParser(formatter_class=argparse.RawDescriptionHelpFormatter,
                                     description=textwrap.dedent("""Python program to collect Dante report tables and 
                                     create overall summary table"""))

    # add arguments
    parser.add_argument('input_path', help='Glob path to directory with input tables.')
    parser.add_argument('output_name', help='Path to output file. Default=\'global_table.tsv\'', nargs='?', default='global_table.tsv')
    parser.add_argument('--sample-regex', '-r', help='Regex for sample name. Default=\'.*\' (all samples)', default='.*')
    parser.add_argument('--table-name', help='Name of input table. Default=\'table.tsv\'', default='table.tsv')
    parser.add_argument('--phasing-table', help='Path to phasing table. Default=\'global_table_phasing.tsv\'', default='global_table_phasing.tsv')

    parsed_args = parser.parse_args()
    parsed_args.output_name = os.path.abspath(parsed_args.output_name)
    parsed_args.phasing_table = os.path.abspath(parsed_args.phasing_table)

    return parsed_args


def create_report(args: argparse.Namespace) -> None:
    # find all motif directories in root directory
    tables = glob.glob(f'{args.input_path}/{args.table_name}')

    # go through all tables and create one global
    global_table = pd.DataFrame()
    for table in tables:
        # if we should skip the sample, skip it:
        sample = table.split('/')[-2]
        if not re.match(args.sample_regex, sample):
            continue

        # read table
        df = pd.read_csv(table, sep='\t')

        # add sample/seq_tech columns
        seq_tech = table.split('/')[-4]
        df['Sample'] = sample
        df['Sequencing Technology'] = seq_tech

        # add to global table
        global_table = pd.concat([global_table, df], ignore_index=True)

    # extract sample number if present
    global_table['Sample ID'] = global_table['Sample'].apply(
        lambda sample_name: int(re.search(r'\d+', sample_name).group(0)) if re.search(r'\d+', sample_name) else None)

    # resort the columns
    first_columns = ['Sequencing Technology', 'Sample', 'Sample ID']
    global_table = global_table[first_columns + [column for column in list(global_table.columns) if column not in first_columns]]
    global_table.sort_values(['Motif', 'Sample ID', 'Sample', 'Sequencing Technology'], ascending=True, inplace=True)

    # split the table
    global_wo_phasing = global_table[~global_table['Repetition index'].str.contains('_', na=True)]
    global_w_phasing = global_table[global_table['Repetition index'].str.contains('_', na=True)]

    # append motifs that have phasing
    motifs_phasing = pd.unique(global_w_phasing['Motif'])
    global_w_phasing = global_table[global_table['Motif'].isin(motifs_phasing)]

    # write the global table to output
    print(f'Writing table with {global_wo_phasing.shape} rows/columns to {args.output_name}')
    global_wo_phasing.to_csv(args.output_name, sep='\t', index=False)
    print(f'Writing table with {global_w_phasing.shape} rows/columns to {args.phasing_table}')
    global_w_phasing.to_csv(args.phasing_table, sep='\t', index=False)


if __name__ == '__main__':
    arguments = load_arguments()
    create_report(arguments)
