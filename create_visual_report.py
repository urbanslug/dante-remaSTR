import argparse
import base64
import enum
import glob
import json
import os
import re
import textwrap
from collections import defaultdict

import pandas as pd
import matplotlib.pyplot as plt

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

    ancestors = parser.add_argument_group('Ancestor graph')

    # add arguments
    parser.add_argument('input_dir', help='Path to directory with Dante reports')
    parser.add_argument('output_dir', help='Path to directory where the output will be stored. Default=<input_dir>/result_files', nargs='?',
                        default=None)
    parser.add_argument('--output-name', '-n', help='Name of the file in the output directory. Default=\'all_visual.html\'',
                        default='all_visual.html')
    parser.add_argument('--sample-regex', '-r', help='Regex for sample name. Default=\'.*\' (all samples)', default='.*')

    ancestors.add_argument('--mother', '-m',
                           help='Name of the sample for mother or 1-based number of the sample lexicographically. Default=1st sample', default=1)
    ancestors.add_argument('--father', '-f',
                           help='Name of the sample for father or 1-based number of the sample lexicographically. Default=2nd sample', default=2)
    ancestors.add_argument('--daughters', '-d',
                           help='Names of samples or their numbers for daughters separated by comma, i.e."1,3" or "wgs_01K,wgs_03K". Default=None',
                           default=None)
    ancestors.add_argument('--sons', '-s',
                           help='Names of samples or their numbers for sons separated by comma, i.e."1,3" or "wgs_01K,wgs_03K". Default=None',
                           default=None)

    args = parser.parse_args()
    if args.output_dir is None:
        args.output_dir = args.input_dir + '/result_files'
    if args.sons is not None:
        args.sons = args.sons.split(',')
    if args.daughters is not None:
        args.daughters = args.daughters.split(',')

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
                        max-width: 600px;
                        word-wrap: break-word;
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
                style = 'style="width: 250px;"' if sample == 'Ancestry' else ''
                cell = f'<img src="data:image/png;base64,{img_data}" {style} alt={str(sample)}>'
            elif str(cell).endswith('.json') and os.path.exists(cell):
                cell = create_plotly(cell, motif_name=(str(index) + str(sample)).replace(' ', ''))
            elif pd.isna(cell):
                cell = ''
            html += f'<td>{cell}</td>\n'
        html += '</tr>\n'
    html += '''
        </table>
        </body>
        </html>
        '''
    return html


def draw_simple_family_tree(mother: str, father: str, daughters: list[str], sons: list[str], filename: str, info: dict = None,
                            wrong_children: list[str] = None) -> None:
    """
    Draw a simple family tree with two parents and a list of children, then save to a file.
    :param mother: Mother name
    :param father: Father name
    :param daughters: List of daughters' names
    :param sons: List of sons' names
    :param filename: Filename to save the image
    :param info: Additional information to display with nodes
    :param wrong_children: List of children that did not pass genotyping check
    """
    if info is None:
        info = {}
    if wrong_children is None:
        wrong_children = []

    fig, ax = plt.subplots()
    ax.axis('off')
    ax.set_xlim(-0.2, 3.1)
    ax.set_ylim(0.6, 2.6)

    # Define coordinates
    m_coord = (0.8, 2)
    f_coord = (2.2, 2)
    chf_coord = (0.5, 1)
    chl_coord = (2.5, 1)

    # Draw template
    ax.plot([m_coord[0], f_coord[0]], [m_coord[1], f_coord[1]], 'k-', linewidth=2)
    ax.plot([(m_coord[0] + f_coord[0]) / 2, (m_coord[0] + f_coord[0]) / 2], [(m_coord[1] + f_coord[1]) / 2, (m_coord[1] + f_coord[1]) / 2 - 0.5],
            'k-', linewidth=2)
    ax.plot([chf_coord[0], chl_coord[0]], [(m_coord[1] + f_coord[1]) / 2 - 0.5, (m_coord[1] + f_coord[1]) / 2 - 0.5], 'k-', linewidth=2)

    # Parent positions
    male_style = 'square,pad=0.4'
    female_style = 'circle,pad=0.3'
    ax.text(*m_coord, f'{mother}\n{info.get(mother, "")}', ha='center', fontsize=20,
            bbox=dict(facecolor='white', edgecolor='black', boxstyle=female_style, linewidth=2))
    ax.text(*f_coord, f'{father}\n{info.get(father, "")}', ha='center', fontsize=20,
            bbox=dict(facecolor='white', edgecolor='black', boxstyle=male_style, linewidth=2))

    # Children
    children = [(d, False) for d in daughters] + [(s, True) for s in sons]
    child_spacing = (chl_coord[0] - chf_coord[0]) / (len(children) - 1)
    for i, (child, male) in enumerate(children):
        # Child coordinates
        x = chf_coord[0] + child_spacing * i
        y = ((len(children) - i - 1) * chf_coord[1] + i * chl_coord[1]) / (len(children) - 1)

        # Draw line to child
        color = 'red' if child in wrong_children else 'black'
        ax.plot([x, x], [(m_coord[1] + f_coord[1]) / 2 - 0.5, y], color=color, linewidth=2)

        # Draw child
        ax.text(x, y, f'{child}\n{info.get(child, "")}', ha='center', fontsize=20,
                bbox=dict(facecolor='white', edgecolor=color, boxstyle=male_style if male else female_style, linewidth=2))

    # Save the figure
    plt.tight_layout()
    plt.savefig(filename, bbox_inches='tight')
    plt.close(fig)


def to_sample_name(sample_id: int | str, samples: list[str]) -> str:
    """
    Search sample name in samples list. If not found try to index by the number.
    :param sample_id: sample number or name
    :param samples: list of samples
    """
    if sample_id in samples:
        return sample_id
    if sample_id.isnumeric():
        return samples[int(sample_id) - 1]
    assert False, sample_id


class ChromEnum(enum.Enum):
    X = 'X'
    Y = 'Y'
    NORM = 'NORM'


def check_genotypes(mother_genotypes: tuple[str, str], father_genotypes: tuple[str, str], child_genotypes: tuple[str, str], son: bool,
                    chrom: ChromEnum = ChromEnum.NORM) -> bool:
    """
    Check if genotypes are valid (children could have been produced by this mother and father genotypes)
    :param mother_genotypes: mother's genotypes
    :param father_genotypes: father's genotypes
    :param child_genotypes: child's genotypes
    :param son: whether this child is son
    :param chrom: which chromosome is this on (normal or X or Y)
    :return: True if genotypes of child are obtainable from mother and father genotypes (or are background or expanded)
    """
    # chrX is inherited from mother for sons and the genotype should be homozygous
    if son and chrom == ChromEnum.X:
        child_genotype = list(set([item for item in child_genotypes if item not in ('B', 'E')]))
        return len(child_genotype) == 0 or (len(child_genotype) == 1 and child_genotype[0] in mother_genotypes)

    # chrY is inherited from father for sons and the genotype should be homozygous
    if son and chrom == ChromEnum.Y:
        child_genotype = list(set([item for item in child_genotypes if item not in ('B', 'E')]))
        return len(child_genotype) == 0 or (len(child_genotype) == 1 and child_genotype[0] in father_genotypes)

    # there should be no calls for chrY motifs in daughters
    if not son and chrom == ChromEnum.Y:
        return child_genotypes[0] == 'B' and child_genotypes[1] == 'B'

    # normal diploid case - one from mother and one from father
    return ((child_genotypes[0] in mother_genotypes + ('B', 'E') and child_genotypes[1] in father_genotypes + ('B', 'E')) or
            (child_genotypes[0] in father_genotypes + ('B', 'E') and child_genotypes[1] in mother_genotypes + ('B', 'E')))


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

    # save the genotyping info for later use
    genotypes = defaultdict(dict)
    nomenclatures = {}

    # fill the table
    table = pd.DataFrame()
    for motif in motifs:
        for sample in samples:
            subtable = pd_tables[sample][pd_tables[sample]['Motif'] == motif]
            nomenclature = subtable['Nomenclature'].iloc[0] if 'Nomenclature' in subtable.columns and len(subtable['Nomenclature']) > 0 else None

            # look if there are more all-call results:
            for i, row in subtable.iterrows():
                number = row['Repetition index']
                repetitions = f'{os.path.realpath(args.input_dir)}/{sample}/{motif}/repetitions_{number}.png'
                if not os.path.exists(repetitions):
                    continue

                parts = number.split('_')
                if len(parts) == 1:
                    number = int(number)
                    allcall = repetitions.replace('repetitions', 'allcall')[:-4] + '.txt'
                    rep_idx2 = None
                    row_name = f'{motif} REP_{number}'
                elif len(parts) == 2:
                    number = int(parts[0])
                    allcall = None
                    rep_idx2 = int(parts[1])
                    row_name = f'{motif} PHASING {number}-{rep_idx2}'
                else:
                    assert False

                # phasing instead of genotyping
                if rep_idx2 is not None:
                    # read phasing:
                    phasing_file = f'{os.path.realpath(args.input_dir)}/{sample}/{motif}/phasing_{number}_{rep_idx2}.txt'
                    phasing_file_contents = load_phasing(phasing_file)
                    (a1, a2), (c, c1, c2) = (('--', '--'), (0.0, 0.0, 0.0)) if phasing_file_contents is None else phasing_file_contents
                    repetitions = repetitions.replace('.png', '.json')
                else:
                    assert os.path.exists(allcall), allcall
                    c, a1, a2, c1, c2, _, _, _, _ = read_all_call(allcall)
                    genotypes[row_name][sample] = (c, a1, a2, c1, c2)

                nomenclatures[row_name] = nomenclature
                table.at[row_name, sample] = repetitions
                sequence = row['Sequence']
                rep_idx = row['Repetition index']
                sequence_highlight = ','.join(
                    [f'<b>{s}</b>' if i + 1 == rep_idx or i + 1 == rep_idx2 else s for i, s in enumerate(sequence.split(','))])
                table.at[row_name + ' result', sample] = (f'alleles: {str(a1):2s} ({float_to_str(c1, percents=True)}) {str(a2):2s} '
                                                          f'({float_to_str(c2, percents=True)}) total {float_to_str(c, percents=True)}<br>'
                                                          f'{sequence_highlight}')

    # add ancestor graphs if needed
    if args.daughters is not None or args.sons is not None:
        # convert numbers to names
        args.mother = to_sample_name(args.mother, samples)
        args.father = to_sample_name(args.father, samples)
        args.daughters = [to_sample_name(d, samples) for d in args.daughters]
        args.sons = [to_sample_name(s, samples) for s in args.sons]

        # go through motifs and fill ancestry
        for row_name in (idx for idx in table.index if 'PHASING' not in idx and 'result' not in idx):
            # define ancestor graph filename
            ancestor_filename = f'{args.input_dir}/{row_name}.png'

            # extract genotyping info
            info = {}
            for sample in args.daughters + args.sons + [args.mother, args.father]:
                c, a1, a2, c1, c2 = genotypes[row_name][sample] if sample in genotypes[row_name] else (0.0, 'B', 'B', 0.0, 0.0)
                info[sample] = f'{a1}/{a2} ({c * 100:.0f}%)'

            # check validity of genotyping info
            mother_genotype = genotypes[row_name][args.mother][1:3] if args.mother in genotypes[row_name] else ('B', 'B')
            father_genotype = genotypes[row_name][args.father][1:3] if args.father in genotypes[row_name] else ('B', 'B')
            nomenclature = nomenclatures[row_name]
            chrom_str = str(nomenclature).split(':')[0].split('.')[0].strip()
            chrom = ChromEnum.X if chrom_str in ['chrX', 'NC_000023'] else (ChromEnum.Y if chrom_str in ['chrY', 'NC_000024'] else ChromEnum.NORM)
            wrong_children = ([child for child in args.daughters if not check_genotypes(mother_genotype, father_genotype,
                                                                                        genotypes[row_name][child][1:3] if child in genotypes[
                                                                                            row_name] else ('B', 'B'), False, chrom)] +
                              [child for child in args.sons if not check_genotypes(mother_genotype, father_genotype,
                                                                                   genotypes[row_name][child][1:3] if child in genotypes[
                                                                                       row_name] else ('B', 'B'), True, chrom)])

            # generate the ancestor graph
            draw_simple_family_tree(args.mother, args.father, args.daughters, args.sons, ancestor_filename, info, wrong_children)

            # append ancestor graph to the table
            table.at[row_name, 'Ancestry'] = ancestor_filename

    # resort the columns
    table = table[(['Ancestry'] + samples) if 'Ancestry' in list(table.columns) else samples]

    # create a html file
    with open(f'{args.output_dir}/{args.output_name}', 'w') as f:
        f.write(df_to_html_table(table))


if __name__ == '__main__':
    args = load_arguments()
    create_report(args)
