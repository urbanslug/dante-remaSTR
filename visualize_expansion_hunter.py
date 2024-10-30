import argparse
from datetime import datetime
import json
import os
import pandas as pd
import textwrap

from src.report.report import plot_histogram_image

html_template = """
<!DOCTYPE html>
<html lang="en">
<head>
    <meta charset="UTF-8">
    <title>Expansion hunter report</title>
    
    <!-- Include Plotly.js -->
    <script src="https://cdn.plot.ly/plotly-latest.min.js"></script>
    <!-- Include jQuery -->
    <script src="https://code.jquery.com/jquery-3.6.0.min.js"></script>
    
    <link rel="stylesheet" type="text/css" href="https://www.w3schools.com/w3css/4/w3.css">

    <script src="https://cdn.datatables.net/1.12.1/js/jquery.dataTables.js"></script>
    <link rel="stylesheet" type="text/css" href="https://cdn.datatables.net/1.12.1/css/jquery.dataTables.css">

    <style>
        body {{
            padding: 20px;
            font-family: Verdana, sans-serif;
            font-size: 15px;
            line-height: 1.5;
        }}

        .tg {{
            border-collapse: collapse;
            border-spacing: 0;
            width: 100%;
        }}

        .tg td {{
            font-weight: normal;
            text-align: center;
            padding: 7px 3px;
            border-style: solid;
            border-width: 1px;
            overflow: hidden;
            word-break: normal;
        }}

        .tg th {{
            font-weight: bold;
            text-align: center;
            padding: 7px 3px;
            border-style: solid;
            border-width: 1px;
            overflow: hidden;
            word-break: normal;
        }}

        table.dataTable {{
            margin: 0;
        }}
    </style>
</head>
<body class="valo v-app">
{body}
</body>
</html>
"""

hist_template = """
<h2 id="{motif_name}">{motif_name}</h2>
alleles: {allele_info}<br>
<div class="" id="hist-{motif_name}"></div>
<script>
    {{
        let hist_data = {motif_reps};

        let updateGraph = () => {{
            hist_data['layout'] = {{...hist_data['layout'], width: (window.innerWidth-50) * 0.6, height: (window.innerWidth-50) * 0.35,
                legend: {{...hist_data['layout']['legend'], x: 0.85, y: 0.85 }}}};
            Plotly.react('hist-{motif_name}', hist_data);
        }};
        updateGraph();
        $(document).ready(function() {{
            window.addEventListener('resize', updateGraph, true);
        }});
    }}
</script>
"""


def load_arguments() -> argparse.Namespace:
    """
    Loads and parses arguments.
    :return: args - parsed arguments
    """
    parser = argparse.ArgumentParser(formatter_class=argparse.RawDescriptionHelpFormatter,
                                     description=textwrap.dedent("""Visualize ExpansionHunter JSON outputs in similar fashion than to Dante."""))

    # add arguments
    parser.add_argument('input_json', help='Path to input JSON file')
    parser.add_argument(
        'output_dir', help='Path to directory where outputs will be stored. Default=<working_dir>', nargs='?', default='eh_viz'
    )

    # parse arguments
    args = parser.parse_args()

    # adjust and create output dir
    if args.output_dir is None:
        args.output_dir = os.getcwd()
    args.output_dir = os.path.abspath(args.output_dir)
    os.makedirs(args.output_dir, exist_ok=True)

    return args


def convert_to_array(text: str) -> list[tuple[int, int]]:
    """
    Convert text representation to array.
    :param text: text representation
    :return: list
    """
    if text == '()':
        return []
    return eval(f'[{text}]')


def search_array(array: list[tuple[int, int]], key: int, up_to: bool = False) -> int:
    """
    Search (rep, count) array for a specific rep and return its count.
    :param array: list[(reps, counts)] - array of repetitions and their counts
    :param key: int - key to search for
    :param up_to: bool - if counting all up to or just equal
    :return: int - count of a specific rep or 0 if not found
    """
    return sum([count for (reps, count) in array if (reps <= key and up_to) or (reps == key and not up_to)])


if __name__ == '__main__':
    # load arguments
    args = load_arguments()

    # print status
    print(f'Visualizing {args.input_json} ... ({datetime.now():%Y-%m-%d %H:%M:%S})')

    # read the json:
    with open(args.input_json) as f:
        data_eh = json.load(f)

    # build table
    table = pd.DataFrame()
    for key, value in data_eh['LocusResults'].items():
        for var_key, var_value in value['Variants'].items():
            motif_clean = var_key.replace('/', '_').replace('+', '_').replace(' ', '_')
            table.at[motif_clean, 'AlleleCount'] = value['AlleleCount']
            table.at[motif_clean, 'Coverage'] = value['Coverage']
            try:
                genotypes = var_value['Genotype'].split('/')
                genotypes_intervals = var_value['GenotypeConfidenceInterval'].split('/')
            except KeyError:
                genotypes = ['0', '0']
                genotypes_intervals = ['0-700', '0-700']
                print('Error in:', args.input_json, key, value)
            table.at[motif_clean, 'Allele 1'] = genotypes[0]
            table.at[motif_clean, 'Allele 2'] = genotypes[1]
            table.at[motif_clean, 'Allele 1 Interval'] = genotypes_intervals[0]
            table.at[motif_clean, 'Allele 2 Interval'] = genotypes_intervals[1]
            table.at[motif_clean, 'RepeatUnit'] = var_value['RepeatUnit']
            table.at[motif_clean, 'ReferenceRegion'] = var_value['ReferenceRegion']
            table.at[motif_clean, 'FigFile'] = f'{args.output_dir}/{motif_clean}'

            # create histograms
            flanking_reads = convert_to_array(var_value['CountsOfFlankingReads'])
            inrepeat_reads = convert_to_array(var_value['CountsOfInrepeatReads'])
            spanning_reads = convert_to_array(var_value['CountsOfSpanningReads'])
            plot_histogram_image(f'{args.output_dir}/{motif_clean}', spanning_reads, flanking_reads, inrepeat_reads)

            # identify allele confirming reads ad count reads
            table.at[motif_clean, 'Reads (flanking)'] = sum(c for (r, c) in flanking_reads)
            table.at[motif_clean, 'Reads (spanning)'] = sum(c for (r, c) in spanning_reads)
            table.at[motif_clean, 'Allele 1 confirming reads (flanking)'] = search_array(flanking_reads, int(genotypes[0]), up_to=True)
            table.at[motif_clean, 'Allele 2 confirming reads (flanking)'] = search_array(flanking_reads, int(genotypes[1]), up_to=True)
            table.at[motif_clean, 'Allele 1 confirming reads (spanning)'] = search_array(spanning_reads, int(genotypes[0]))
            table.at[motif_clean, 'Allele 2 confirming reads (spanning)'] = search_array(spanning_reads, int(genotypes[1]))
            table.at[motif_clean, 'Reads (inrepeat)'] = sum(c for (r, c) in inrepeat_reads)

    # plotting figures
    print(f'Plotting figures in {args.output_dir} ... ({datetime.now():%Y-%m-%d %H:%M:%S})')

    # plot table and all the result histograms
    body = table.drop('FigFile', axis=1).to_html(classes='tg')
    for i, row in table.iterrows():
        motif_reps = open(f'{row["FigFile"]}.json').read()
        allele_info = f'{row["Allele 1"]} / {row["Allele 2"]} ({row["Allele 1 Interval"]} / {row["Allele 2 Interval"]})'
        histogram = hist_template.format(motif_name=i, motif_reps=motif_reps, allele_info=allele_info)
        body += '\n' + histogram

    # write into a html file
    with open(f'{args.output_dir}/report.html', 'wt') as f:
        print(html_template.format(body=body), file=f)

    # write into a tsv table
    int_columns = ['AlleleCount', 'Reads (flanking)', 'Reads (spanning)', 'Reads (inrepeat)', 'Allele 1 confirming reads (flanking)',
                   'Allele 2 confirming reads (flanking)', 'Allele 1 confirming reads (spanning)', 'Allele 2 confirming reads (spanning)']
    table.drop('FigFile', axis=1).astype({k: int for k in int_columns}).to_csv(f'{args.output_dir}/results.tsv', sep='\t')
