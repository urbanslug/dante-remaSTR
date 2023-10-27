import argparse
import base64
import glob
import os
import textwrap

import pandas as pd

from src.report.report import read_all_call


def load_arguments() -> argparse.Namespace:
    """
    Loads and parses arguments.
    :return: args - parsed arguments
    """
    parser = argparse.ArgumentParser(formatter_class=argparse.RawDescriptionHelpFormatter,
                                     description=textwrap.dedent("""Python program to collect Dante histogram files and 
                                     create 'full' report, to make the comparison of multiple samples and motifs easier"""))

    # add arguments
    parser.add_argument('input_dir', help='Path to directory with Dante reports')
    parser.add_argument('output_dir', help='Path to directory where the output will be stored. Default=<input_dir>/result_files', nargs='?',
                        default=None)

    args = parser.parse_args()
    if args.output_dir is None:
        args.output_dir = args.input_dir + '/result_files'

    return args


def encode_image(image_path):
    """Encode image to base64"""
    with open(image_path, 'rb') as img_file:
        return base64.b64encode(img_file.read()).decode('utf-8')


def df_to_html_table(df):
    """converts a DataFrame to an HTML table without borders and with full-width images"""
    col_count = len(df.columns)
    img_width = 100 / col_count
    html = '<table style="border-collapse: collapse; width: 200%;">'
    # header
    html += '<tr>'
    html += f'<th style="border: none; vertical-align: middle; text-align: center;">Motif</th>'
    for col in df.columns:
        html += f'<th style="border: none; vertical-align: middle; text-align: center;">{col}</th>'
    html += '</tr>'
    # rows
    for i, row in df.iterrows():
        html += '<tr>'
        html += f'<td style="border: none;"><strong>{i}</strong></td>'  # print row name in bold
        for cell in row:
            if str(cell).endswith(('.png', '.jpg', '.jpeg')) and os.path.exists(cell):
                img_data = encode_image(cell)
                cell = f'<img src="data:image/png;base64,{img_data}" style="width: 100%; height: auto;" />'
            html += f'<td style="border: none; vertical-align: middle; text-align: center;">{cell}</td>'
        html += '</tr>'
    html += '</table>'
    return html


def create_report(args: argparse.Namespace):
    # create output directory
    os.makedirs(args.output_dir, exist_ok=True)

    # find all motif directories in root directory
    possible_paths = glob.glob(f'{args.input_dir}/*/*/allcall_1.txt')
    samples = sorted(set([path.split('/')[-3] for path in possible_paths]))
    motifs = sorted(set([path.split('/')[-2] for path in possible_paths]))

    # read all result tables
    pd_tables = {sample: pd.read_csv(f'{args.input_dir}/{sample}/table.tsv', sep='\t') for sample in samples}

    # fill the table
    table = pd.DataFrame()
    for motif in motifs:
        for sample in samples:
            # look if there are more allcall results:
            for allcall in sorted(glob.glob(f'{os.path.realpath(args.input_dir)}/{sample}/{motif}/allcall_*.txt')):
                repetitions = allcall.replace('allcall', 'repetitions')[:-4] + '.png'
                number = allcall[:-4].split('_')[-1]
                table.at[motif + f' REP.{number}', sample] = repetitions
                c, a1, a2, c1, c2, _, _, _, _ = read_all_call(allcall)
                table_row = pd_tables[sample][pd_tables[sample]['Motif'] == motif].iloc[int(number) - 1]
                sequence = table_row['Sequence']
                rep_idx = table_row['Repetition index']
                sequence_highlight = ','.join([f'<b>{s}</b>' if i == rep_idx - 1 else s for i, s in enumerate(sequence.split(','))])
                table.at[motif + f' REP.{number} result', sample] = (f'alleles: {str(a1):2s} ({c1 * 100: 5.1f}%) {str(a2):2s} '
                                                                     f'({c2 * 100: 5.1f}%) total {c * 100: 5.1f}%</br>'
                                                                     f'{sequence_highlight}')

    # go through them and create a html file
    with open(f'{args.output_dir}/all_visual.html', 'w') as f:
        f.write(df_to_html_table(table))


if __name__ == '__main__':
    args = load_arguments()
    create_report(args)
