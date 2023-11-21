import argparse
import base64
import glob
import json
import os
import re
import textwrap

import pandas as pd

from src.inference import load_phasing
from src.report.html_templates import float_to_str
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
    parser.add_argument('--output-name', '-n', help='Name of the file in the output directory. Default=\'all_visual.html\'',
                        default='all_visual.html')
    parser.add_argument('--sample-regex', '-s', help='Regex for sample name. Default=\'.*\' (all samples)', default='.*')

    args = parser.parse_args()
    if args.output_dir is None:
        args.output_dir = args.input_dir + '/result_files'

    return args


def encode_image(image_path: str) -> str:
    """Encode image to base64"""
    with open(image_path, 'rb') as img_file:
        return base64.b64encode(img_file.read()).decode('utf-8')


def create_plotly(json_path: str, motif_name: str) -> str:
    """Create plotly image"""

    # define template
    histogram_template = """<div id="hist-{motif_name}" style="width: 600px; height: 500px;"></div>
            <script>
                Plotly.newPlot('hist-{motif_name}', {motif_reps});
            </script>"""

    # read the plotly file
    with open(json_path, 'r') as json_file:
        # json_data = json_file.read()
        json_data = json.load(json_file)
        json_data['layout']['width'] = '600'
        json_data['layout']['height'] = '500'

    # fill template and return
    return histogram_template.format(motif_name=motif_name, motif_reps=json.dumps(json_data))


def df_to_html_table(df: pd.DataFrame) -> str:
    """Converts a DataFrame to an HTML table without borders and with full-width images"""
    html = '''
        <!DOCTYPE html>
        <html lang="en">
            <head>
                <meta charset="UTF-8">
                <title>Table Export</title>
                <script src="https://ajax.googleapis.com/ajax/libs/jquery/3.6.1/jquery.min.js"></script>
                <script src="https://cdn.plot.ly/plotly-2.14.0.min.js"></script>
                <style>
                    table, th, td {
                        border-collapse: collapse;
                        border: none;
                        vertical-align: middle;
                        text-align: center;
                    }
                    th, td {
                        padding: 8px;
                    }
                    img {
                        width: 600px; 
                        height: auto;
                    }
                    </style>
            </head>
        <body>
            <table>
        '''
    # header
    html += '<tr>'
    html += f'<th>Motif</th>'
    for col in df.columns:
        html += f'<th>{col}</th>'
    html += '</tr>\n'
    # rows
    for index, row in df.iterrows():
        html += '<tr>'
        html += f'<td><strong>{index}</strong></td>\n'  # print row name in bold
        for sample, cell in row.items():
            if str(cell).endswith(('.png', '.jpg', '.jpeg')) and os.path.exists(cell):
                img_data = encode_image(cell)
                cell = f'<img src="data:image/png;base64,{img_data}">'
            if str(cell).endswith('.json') and os.path.exists(cell):
                cell = create_plotly(cell, motif_name=(str(index) + str(sample)).replace(' ', ''))
            html += f'<td>{cell}</td>\n'
        html += '</tr>\n'
    html += '''
        </table>
        </body>
        </html>
        '''
    return html


def create_report(args: argparse.Namespace) -> None:
    # create output directory
    os.makedirs(args.output_dir, exist_ok=True)

    # find all motif directories in root directory
    possible_paths = glob.glob(f'{args.input_dir}/*/*/allcall_1.txt')
    samples = set([path.split('/')[-3] for path in possible_paths])
    samples = sorted([sample for sample in samples if re.match(args.sample_regex, sample)])
    motifs = sorted(set([path.split('/')[-2] for path in possible_paths]))

    # read all result tables
    pd_tables = {sample: pd.read_csv(f'{args.input_dir}/{sample}/table.tsv', sep='\t') for sample in samples}

    # fill the table
    table = pd.DataFrame()
    for motif in motifs:
        for sample in samples:
            subtable = pd_tables[sample][pd_tables[sample]['Motif'] == motif]
            single_reps = len(subtable)
            # look if there are more all-call results:
            for repetitions in sorted(glob.glob(f'{os.path.realpath(args.input_dir)}/{sample}/{motif}/repetitions_*.png')):
                parts = os.path.basename(repetitions[:-4]).split('_')
                if len(parts) == 2:
                    number = int(parts[1])
                    allcall = repetitions.replace('repetitions', 'allcall')[:-4] + '.txt'
                    rep_idx2 = None
                    row_name = f'{motif} REP.{number}'
                elif len(parts) == 3:
                    number = int(parts[1])
                    allcall = None
                    rep_idx2 = int(parts[2])
                    row_name = f'{motif} PHASING {number}-{number + 1}'
                else:
                    assert False

                # phasing instead of genotyping TODO add numeric phasing
                if number > single_reps:
                    number -= single_reps
                    row_name = f'{motif} PHASING {number}-{number + 1}'
                    rep_idx2 = subtable.iloc[number]['Repetition index']
                    c, a1, a2, c1, c2 = 0.0, '--', '--', 0.0, 0.0
                    repetitions = repetitions.replace('.png', '.json')
                elif rep_idx2 is not None:
                    # read phasing:
                    phasing_file = f'{os.path.realpath(args.input_dir)}/{sample}/{motif}/phasing_{number}_{rep_idx2}.txt'
                    phasing_file_contents = load_phasing(phasing_file)
                    (a1, a2), (c, c1, c2) = (('--', '--'), (0.0, 0.0, 0.0)) if phasing_file_contents is None else phasing_file_contents
                    repetitions = repetitions.replace('.png', '.json')
                else:
                    assert os.path.exists(allcall), allcall
                    c, a1, a2, c1, c2, _, _, _, _ = read_all_call(allcall)

                table.at[row_name, sample] = repetitions
                table_row = subtable.iloc[number - 1]
                sequence = table_row['Sequence']
                rep_idx = table_row['Repetition index']
                sequence_highlight = ','.join(
                    [f'<b>{s}</b>' if i + 1 == rep_idx or i + 1 == rep_idx2 else s for i, s in enumerate(sequence.split(','))])
                table.at[row_name + ' result', sample] = (f'alleles: {str(a1):2s} ({float_to_str(c1, percents=True)}) {str(a2):2s} '
                                                          f'({float_to_str(c2, percents=True)}) total {float_to_str(c, percents=True)}<br />'
                                                          f'{sequence_highlight}')
    # resort the columns
    table = table[samples]

    # go through them and create a html file
    with open(f'{args.output_dir}/{args.output_name}', 'w') as f:
        f.write(df_to_html_table(table))


if __name__ == '__main__':
    args = load_arguments()
    create_report(args)
