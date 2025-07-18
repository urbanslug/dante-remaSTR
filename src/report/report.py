import enum
import os
import re
import shutil
import typing
import gzip
from collections import Counter

import matplotlib
import plotly.graph_objects as go
import numpy as np
import pandas as pd

matplotlib.use('Agg')
import matplotlib.pyplot as plt

from src import annotation
from src import report
from src.annotation.Motif import Motif
from src.postfilter import PostFilter

# max repetitions on the graph
MAX_REPETITIONS = 40


class ChromEnum(enum.Enum):
    X = 'X'
    Y = 'Y'
    NORM = 'NORM'


def chrom_from_string(chrom_str: str) -> ChromEnum:
    """
    Converts a string to a ChromEnum object.
    :param chrom_str: str - the string to convert to a ChromEnum object
    :return ChromEnum - enum object representing the chromosome
    """
    return ChromEnum.X if chrom_str in ['chrX', 'NC_000023'] else (ChromEnum.Y if chrom_str in ['chrY', 'NC_000024'] else ChromEnum.NORM)


def write_annotations(out_file: str, annotations: list[annotation.Annotation], zip_it: bool = True) -> None:
    """
    Stores annotations in alignment format into output text file
    :param out_file: Alignment file
    :param annotations: Annotated reads
    :param zip_it: bool - whether to gzip the resulting file
    """
    if zip_it and not out_file.endswith('.gz'):
        out_file += '.gz'
    with gzip.open(out_file, 'wt') if zip_it else open(out_file, 'w') as fw:
        for annot in annotations:
            fw.write(str(annot) + '\n')


def write_annotation_pairs(out_file: str, annotation_pairs: list[annotation.AnnotationPair], zip_it: bool = True) -> None:
    """
    Stores annotations in alignment format into output text file
    :param out_file: Alignment file
    :param annotation_pairs: Annotated pairs of reads
    :param zip_it: bool - whether to gzip the resulting file
    """
    if zip_it and not out_file.endswith('.gz'):
        out_file += '.gz'
    with gzip.open(out_file, 'wt') if zip_it else open(out_file, 'w') as fw:
        for ap in annotation_pairs:
            write_left = str(ap.ann1) if ap.ann1 is not None else 'Left None\n'
            fw.write(write_left)
            write_right = str(ap.ann2) if ap.ann2 is not None else 'Right None\n'
            fw.write(write_right + '\n')


def write_alignment(out_file: str, annotations: list[annotation.Annotation], index_rep: int, allele: int = None, index_rep2: int = None,
                    allele2: int = None, zip_it: bool = True, cutoff_after: int = None, right_align: bool = False) -> None:
    # TODO this needs complete rework
    """
    Creates a multi-alignment of all annotations into output text file
    :param out_file: str - alignment filename
    :param annotations: list(Annotation) - annotated reads
    :param index_rep: int - index of repetition module of a motif
    :param allele: int/None - which allele to print only, if None print all of them
    :param index_rep2: int - index of second repetition module of a motif
    :param allele2: int/None - which allele2 to print only, if None print all of them
    :param zip_it: bool - whether to gzip the resulting file
    :param right_align: bool - whether we deal with right alignment file
    :param cutoff_after: int - how many bases to keep outside the annotated motif (None for keep all)
    """
    # select annotations
    if allele is not None:
        if allele2 is not None and index_rep2 is not None:
            annotations = [a for a in annotations if a.module_repetitions[index_rep] == allele and a.module_repetitions[index_rep2] == allele2]
        else:
            annotations = [a for a in annotations if a.module_repetitions[index_rep] == allele]

    # apply cutoff
    if cutoff_after is not None:
        annotations = [a.get_shortened_annotation(cutoff_after) for a in annotations]

    # setup alignments
    alignments = [''] * len(annotations)  # alignment strings
    align_inds = np.zeros(len(annotations), dtype=int)  # indices of annotations that were processed
    states = []  # has numbers of states in the final multiple alignment

    while True:
        # get minimal state:
        min_comp = (True, 'Z', -1)
        total_done = 0
        for i, (annot, ai) in enumerate(zip(annotations, align_inds)):
            if ai >= len(annot.states):
                total_done += 1
                continue
            state = annot.states[ai]
            comparator = (state != 'I', state, i)
            min_comp = min(comparator, min_comp)

        # if we have done every state, end:
        if total_done >= len(alignments):
            break

        states.append(min_comp[1])

        # now print all states, that are minimal:
        for i, (annot, ai) in enumerate(zip(annotations, align_inds)):
            if ai >= len(annot.states):
                alignments[i] += '_'  # put ends of each alignment to be of same length
                continue
            if annot.states[ai] == min_comp[1]:
                alignments[i] += annot.read_seq[ai]
                align_inds[i] += 1
            else:
                alignments[i] += '_'

    # sort according to motif count:
    left_flank = np.array([-(ann.module_bases[0] + ann.left_flank_len) for ann in annotations])
    left_flank_exist = np.array([-(ann.module_repetitions[0]) if not right_align else 0 for ann in annotations])
    if index_rep2 is not None:
        # sorting first with 1st allele then with second
        reps1 = np.array([-ann.module_bases[index_rep] for ann in annotations])
        reps2 = np.array([-ann.module_bases[index_rep2] for ann in annotations])
        sort_inds = np.lexsort((left_flank, reps2, reps1, left_flank_exist))  # sort by existence of left flank, first allele, second, left flank len.
    else:
        reps = np.array([-ann.module_bases[index_rep] for ann in annotations])
        sort_inds = np.lexsort((left_flank, reps, left_flank_exist))
    annotations = np.array(annotations)[sort_inds]
    alignments = list(np.array(alignments)[sort_inds])

    def move_right(alignment: str, start: int, end: int) -> str:
        """
        Shift first part of the alignment to the right.
        :param alignment: str - alignment of the read
        :param start: int - start idx for shift
        :param end: int - one after end idx for a shift
        :return: str - alignment, where first part is shifted to right
        """
        align_part = alignment[start:end]
        # find last empty:
        idx = 0
        for idx in reversed(range(len(align_part))):
            if align_part[idx] != '_':
                break
        idx += 1

        # return shifted alignment
        return alignment[:start] + ('_' * (len(align_part) - idx)) + align_part[:idx] + alignment[end:]

    def get_range(symbol: str = '0') -> (int, int):
        """
        Get range of a module
        :param symbol: str - symbol for state to get range for
        :return: int, int - start and (one after) end range coordinates
        """
        try:
            first_idx = states.index(symbol)
            last_idx = len(states) - states[-1::-1].index(symbol) - 1
            return first_idx, last_idx + 1
        except ValueError:
            return -1, -1

    def get_left_flank() -> int:
        """
        Get range of a left flank.
        :return: int - (one after) nd of the left flank before module '0'
        """
        for i, state in enumerate(states):
            if state != '-':
                return i
        return -1

    # for every alignment, shift the left flank right
    end = get_left_flank()
    if end != -1:
        for i in range(len(alignments)):
            alignments[i] = move_right(alignments[i], 0, end)

    # for every alignment, shift the first module right
    start0, end0 = get_range('0')
    if start0 != -1:
        for i in range(len(alignments)):
            alignments[i] = move_right(alignments[i], start0, end0)

    # in addition, those that have only '_' in state '0' (missing left flank), shift right also '1' state
    start1, end1 = get_range('1')
    first_zero_idx = len(alignments)
    for i in range(len(alignments)):
        if start1 != -1 and (start0 == -1 or alignments[i][start0:end0].count('_') == end0 - start0 or right_align):
            first_zero_idx = min(first_zero_idx, i)
            alignments[i] = move_right(alignments[i], start1, end1)

    # add empty line if we have some alignments without left flank
    annot_names = [annot.read_id for annot in annotations]
    if first_zero_idx != len(alignments) and not right_align:
        alignments = alignments[:first_zero_idx] + ['_' * len(alignments[0])] + alignments[first_zero_idx:]
        annot_names = annot_names[:first_zero_idx] + ['empty_line'] + annot_names[first_zero_idx:]

    # print to file
    if zip_it and not out_file.endswith('.gz'):
        out_file += '.gz'
    with gzip.open(out_file, 'wt') if zip_it else open(out_file, 'w') as fw:
        # print alignments
        for annot_name, align in zip(annot_names, alignments):
            print(f'>{annot_name}', file=fw)
            print(align, file=fw)


def sorted_repetitions(annotations: list[annotation.Annotation]) -> list[tuple[tuple[int, ...], int]]:
    """
    Aggregate same repetition counts for annotations and sort them according to quantity of repetitions of each module
    :param annotations: Annotated reads
    :return: list of (repetitions, count), sorted by repetitions
    """
    count_dict = Counter(tuple(annot.module_repetitions) for annot in annotations)
    return sorted(count_dict.items(), key=lambda k: k[0])


def write_histogram(out_file: str, annotations: list[annotation.Annotation], ) -> None:
    """
    Stores quantity of different combinations of module repetitions into text file
    :param out_file: str - output file for repetitions
    :param annotations: Annotated reads
    """
    # setup
    sorted_reps = sorted_repetitions(annotations)

    # write repetitions.txt
    with open(out_file, 'w') as fw:
        for repetitions, counts in sorted_reps:
            rep_code = '\t'.join(map(str, repetitions))
            fw.write(f'{counts}\t{rep_code}\n')


def write_profile(profile_file: str, annotations: list[annotation.Annotation], index_rep: int) -> None:
    """
    Stores quantity of different combinations of module repetitions into text file
    :param profile_file: str - output file for repetitions
    :param annotations: Annotated reads
    :param index_rep: int - index of the first repetition
    """
    # setup
    sorted_reps = sorted_repetitions(annotations)

    # write profile
    length = max([0] + [x[0][index_rep] for x in sorted_reps])
    profile = np.zeros(length + 1, dtype=int)

    for repetitions, counts in sorted_reps:
        profile[repetitions[index_rep]] += counts

    with open(profile_file, 'w') as f:
        f.write('\t'.join(map(str, profile)))


def write_histogram_nomenclature(out_file: str, annotations: list[annotation.Annotation],
                                 index_rep: int = None, index_rep2: int = None) -> None:
    """
    Stores quantity of different nomenclature strings into text file
    :param out_file: str - output file for repetitions
    :param annotations: Annotated reads
    :param index_rep: int - index of the first repetition (None if include all)
    :param index_rep2: int - index of the second repetition (None if include all)
    """
    # count nomenclature strings:
    count_dict = Counter(annot.get_nomenclature(index_rep, index_rep2, False) for annot in annotations)
    count_dict = sorted(count_dict.items(), key=lambda k: (-k[1], k[0]))

    # write nomenclatures to file
    with open(out_file, 'w') as fw:
        for nomenclature, count in count_dict:
            fw.write(f'{count}\t{nomenclature}\n')


def write_histogram_image2d(out_prefix: str, deduplicated: list[annotation.Annotation],
                            index_rep: int, index_rep2: int, seq: str, seq2: str) -> None:
    """
    Stores quantity of different combinations of module repetitions, generates separate graph image for each module
    :param out_prefix: Output file prefix
    :param deduplicated: list[Annotation] - read pairs
    :param index_rep: int - index of repetition module of a motif
    :param index_rep2: int - index of the second repetition module of a motif
    :param seq: str - module of the repetition
    :param seq2: str - 2nd module of the repetition
    """
    if deduplicated is None or len(deduplicated) == 0:
        return

    dedup_reps = [(x.get_str_repetitions(index_rep), x.get_str_repetitions(index_rep2)) for x in deduplicated
                  if x.get_str_repetitions(index_rep) is not None and x.get_str_repetitions(index_rep2) is not None]

    if len(dedup_reps) == 0:
        return

    # assign maximals
    xm = max([r for (_, r), _ in dedup_reps])
    ym = max([r for _, (_, r) in dedup_reps])
    max_ticks = max(ym, xm) + 2
    xm = max(MAX_REPETITIONS, xm)
    ym = max(MAX_REPETITIONS, ym)

    # create data containers
    data = np.zeros((xm + 1, ym + 1), dtype=int)
    data_primer = np.zeros((xm + 1, ym + 1), dtype=int)
    for (c1, r1), (c2, r2) in dedup_reps:
        if c1 and c2:
            data[r1, r2] += 1
        if c1 and not c2:
            data_primer[r1, r2:] += 1
        if not c1 and c2:
            data_primer[r1:, r2] += 1

    # create colormaps:
    cmap_blue = matplotlib.cm.get_cmap('Blues')
    cmap_blue = cmap_blue(np.arange(int(cmap_blue.N * 0.15), int(cmap_blue.N * 0.8)))  # start from light blue to deep blue
    cmap_blue[0, -1] = 0.0  # Set alpha on the lowest element only
    cmap_grey = matplotlib.cm.get_cmap('Greys')
    cmap_grey = cmap_grey(np.arange(int(cmap_grey.N * 0.15), int(cmap_grey.N * 0.6)))  # start from light grey to deep grey
    cmap_grey[0, -1] = 0.0  # Set alpha on the lowest element only
    cmap_grey_plotly = [(i, f'rgba({c[0]}, {c[1]}, {c[2]}, {c[3]})') for i, c in [(0.0, cmap_grey[0]), (0.01, cmap_grey[1]), (1.0, cmap_grey[-1])]]
    cmap_blue_plotly = [(i, f'rgba({c[0]}, {c[1]}, {c[2]}, {c[3]})') for i, c in [(0.0, cmap_blue[0]), (0.01, cmap_blue[1]), (1.0, cmap_blue[-1])]]

    # plot pcolor
    plt.figure(figsize=(12, 8))
    img2 = plt.pcolor(data_primer[:max_ticks, :max_ticks], cmap=matplotlib.colors.ListedColormap(cmap_grey), alpha=0.4,
                      edgecolor=(1.0, 1.0, 1.0, 0.0), lw=0, vmin=np.min(data_primer),
                      vmax=np.max(data_primer) + 0.01)
    img1 = plt.pcolor(data[:max_ticks, :max_ticks], cmap=matplotlib.colors.ListedColormap(cmap_blue), vmin=np.min(data), vmax=np.max(data) + 0.01)
    plt.xticks()
    plt.ylabel('STR %d [%s]' % (index_rep + 1, seq.split('-')[-1]))
    plt.xlabel('STR %d [%s]' % (index_rep2 + 1, seq2.split('-')[-1]))
    plt.colorbar(img1)
    plt.colorbar(img2)

    # setup ticks
    start_ticks = 5
    step_ticks = 5
    plt.xticks(np.array(range(start_ticks, max_ticks + 1, step_ticks)) + 0.5, range(start_ticks, max_ticks + 1, step_ticks))
    plt.yticks(np.array(range(start_ticks, max_ticks + 1, step_ticks)) + 0.5, range(start_ticks, max_ticks + 1, step_ticks))

    # output it
    plt.savefig(out_prefix + '.pdf')
    plt.savefig(out_prefix + '.png')
    plt.close()

    # ----- PLOTLY HISTOGRAM -----
    def parse_labels(num, num_primer):
        if num == 0 and num_primer == 0:
            return ''
        elif num == 0 and num_primer != 0:
            return '0/%s' % str(num_primer)
        elif num != 0 and num_primer == 0:
            return '%s/0' % str(num)
        else:
            return '%s/%s' % (str(num), str(num_primer))

    str1 = 'STR %d [%s]' % (index_rep + 1, seq.split('-')[-1])
    str2 = 'STR %d [%s]' % (index_rep2 + 1, seq2.split('-')[-1])

    text = [[parse_labels(data[i, j], data_primer[i, j]) for j in range(data.shape[1])] for i in range(data.shape[0])]

    fig = go.Figure()
    if np.sum(data_primer[:max_ticks, :max_ticks]) > 0:
        fig.add_trace(go.Heatmap(z=data_primer[:max_ticks, :max_ticks], name='Repetitions heatmap',
                                 showscale=True, colorbar_x=1.3, colorbar_title='Partial reads', colorscale=cmap_grey_plotly))
    fig.add_trace(go.Heatmap(z=data[:max_ticks, :max_ticks], text=text, name='Repetitions heatmap',
                             showscale=True, colorbar_title='Full reads', colorscale=cmap_blue_plotly))

    fig.update_traces(texttemplate='%{text}', textfont_size=7,
                      hovertemplate='<b>{name1}:\t%{y}<br>{name2}:\t%{x}</b><br>Full / Partial:\t%{text}'.
                      format(name1=str1, y='{y}', name2=str2, x='{x}', text='{text}'))
    fig.update_layout(width=800, height=600, template='simple_white')
    fig.update_yaxes(title_text=str1)
    fig.update_xaxes(title_text=str2)

    with open(out_prefix + '.json', 'w') as f:
        f.write(fig.to_json())

    # fig.write_image(out_prefix + '_plotly.pdf')


def write_histogram_image(out_prefix: str, annotations: list[annotation.Annotation],
                          filt_annot: list[annotation.Annotation], index_rep: int) -> None:
    """
    Stores quantity of different combinations of module repetitions, generates separate graph image for each module
    :param out_prefix: Output file prefix
    :param annotations: Annotated reads.
    :param filt_annot: Annotated reads (filtered)
    :param index_rep: int - index of repetition module of a motif
    """
    if len(annotations) == 0 and len(filt_annot) == 0:
        return

    repetitions = sorted_repetitions(annotations)
    repetitions_filt = sorted_repetitions(filt_annot)

    plot_histogram_image(out_prefix, [(r[index_rep], c) for r, c in repetitions], [(r[index_rep], c) for r, c in repetitions_filt])


def plot_histogram_image(out_prefix: str, spanning_counts: list[tuple[int, int]], filtered_counts: list[tuple[int, int]],
                         inread_counts: list[tuple[int, int]] = None) -> None:
    """
    Generates separate graph image for each module
    :param out_prefix: Output file prefix
    :param spanning_counts: list of spanning counts (repetition, number of reads)
    :param filtered_counts: list of flanking counts (repetition, number of reads)
    :param inread_counts: list of inread counts (repetition, number of reads)
    """
    # empty inread array
    if inread_counts is None:
        inread_counts = []

    # adjust variables
    width = 0.9
    plt.figure(figsize=(20, 8))
    xm = max([r for r, c in spanning_counts] + [r for r, c in filtered_counts] + [r for r, c in inread_counts] + [MAX_REPETITIONS])
    dist = [0] * (xm + 1)

    # set data
    for r, c in spanning_counts:
        dist[r] += c
    dist_filt = dist.copy()
    for r, c in filtered_counts:
        dist_filt[r] += c
    dist_inread = dist_filt.copy()
    for r, c in inread_counts:
        dist_inread[r] += c

    # create barplots
    rects_inread = [None] * len(dist)
    if len(inread_counts) > 0:
        rects_inread = plt.bar(np.arange(xm + 1), dist_inread, width, color='orange', alpha=0.4)
    rects_filt = plt.bar(np.arange(xm + 1), dist_filt, width, color='lightgrey')
    rects = plt.bar(np.arange(xm + 1), dist, width)
    plt.xticks(np.arange(1, xm + 1))
    plt.ylabel('Counts')
    plt.xlabel('STR repetitions')
    _, max_y = plt.ylim()
    plt.xlim((0, xm + 1))

    # label numbers
    for rect, rect_filt, rect_inread in zip(rects, rects_filt, rects_inread):
        if rect.get_height() > 0:
            plt.text(rect.get_x() + rect.get_width() / 2., rect.get_height() + max_y / 100.0, '%d' % int(rect.get_height()), ha='center', va='bottom')
        if rect_filt.get_height() != rect.get_height():
            plt.text(rect_filt.get_x() + rect_filt.get_width() / 2., rect_filt.get_height() + max_y / 100.0,
                     '%d' % int(rect_filt.get_height() - rect.get_height()), ha='center', va='bottom', color='grey')
        if rect_inread is not None and rect_inread.get_height() != rect_filt.get_height():
            plt.text(rect_inread.get_x() + rect_inread.get_width() / 2., rect_inread.get_height() + max_y / 100.0,
                     '%d' % int(rect_inread.get_height() - rect_filt.get_height()), ha='center', va='bottom', color='orange')

    # output it
    plt.savefig(out_prefix + '.pdf')
    plt.savefig(out_prefix + '.png')
    plt.close()

    # ----- PLOTLY HISTOGRAM -----
    dist_text = ['' if d == 0 else str(d) for d in dist]
    dist_filt_text = ['' if df - d == 0 else str(df - d) for df, d in zip(dist_filt, dist)]
    dist_inread_text = ['' if di - df == 0 else str(di - df) for di, df in zip(dist_inread, dist_filt)]

    fig = go.Figure()
    if len(inread_counts) > 0:
        fig.add_bar(y=dist_inread, text=dist_inread_text, name='Inread reads', marker_color='#FF6600', textfont_color='#FF6600')
    fig.add_bar(y=dist_filt, text=dist_filt_text, name='Partial reads', marker_color='#CCCCCC', textfont_color='#CCCCCC')
    fig.add_bar(y=dist, text=dist_text, name='Full reads', marker_color='#636EFA', textfont_color='#636EFA')

    fig.update_traces(textposition='outside', texttemplate='%{text}', hovertemplate='%{text}', textfont_size=7)
    fig.update_layout(width=800, height=450,
                      title='Histogram of repetitions',
                      hovermode='x',
                      yaxis_fixedrange=True,
                      template='simple_white',
                      barmode='overlay')
    fig.update_yaxes(title_text='Read counts')
    fig.update_xaxes(title_text='STR repetitions', tickmode='array',
                     tickvals=list(range(5, len(dist), 5)),
                     ticktext=list(range(5, len(dist), 5)))

    with open(out_prefix + '.json', 'w') as f:
        f.write(fig.to_json())

    # fig.write_image(out_prefix + '_plotly.pdf')


def write_all(quality_annotations: list[annotation.Annotation], filt_primer: list[annotation.Annotation],
              filtered_annotations: list[annotation.Annotation], motif_dir: str, motif_class: Motif, module_number: int,
              second_module_number: int | None = None, zip_it: bool = True, cutoff_alignments: int | None = None) -> None:
    """
    Write all output files: quality annotations, one-primer annotations, filtered annotations, statistics, repetitions + images.
    :param quality_annotations: list(Annotation) - list of blue annotations
    :param filt_primer: list(Annotation) - list of grey annotations
    :param filtered_annotations: list(Annotation) - list of filtered out annotation
    :param motif_dir: str - path to motif directory
    :param motif_class: Motif - motif class
    :param module_number: int - index of first studied repetition in modules
    :param second_module_number: int - index of second studied repetition in modules (optional)
    :param zip_it: bool - whether to gzip the resulting file
    :param cutoff_alignments: int - how many bases to keep beyond annotated modules in alignments
    :return: None
    """
    # create dir if not exists
    os.makedirs(motif_dir, exist_ok=True)

    # create suffix for files:
    suffix = str(module_number) if second_module_number is None else f'{module_number}_{second_module_number}'

    # write output files
    write_annotations(f'{motif_dir}/annotations_{suffix}.txt', quality_annotations, zip_it=zip_it)
    if len(quality_annotations) and max([len(a.read_seq) for a in quality_annotations]) > 200:
        write_annotations(f'{motif_dir}/annotations_{suffix}_short.txt', [a.get_shortened_annotation(cutoff_alignments) for a in quality_annotations],
                          zip_it=zip_it)
    write_annotations(f'{motif_dir}/filtered_{suffix}.txt', filtered_annotations, zip_it=zip_it)
    if len(filtered_annotations) and max([len(a.read_seq) for a in filtered_annotations]) > 200:
        write_annotations(f'{motif_dir}/filtered_{suffix}_short.txt', [a.get_shortened_annotation(cutoff_alignments) for a in filtered_annotations],
                          zip_it=zip_it)
    write_annotations(f'{motif_dir}/filtered_primer_{suffix}.txt', filt_primer, zip_it=zip_it)
    if len(filt_primer) and max([len(a.read_seq) for a in filt_primer]) > 200:
        write_annotations(f'{motif_dir}/filtered_primer_{suffix}_short.txt', [a.get_shortened_annotation(cutoff_alignments) for a in filt_primer],
                          zip_it=zip_it)
    write_alignment(f'{motif_dir}/alignment_{suffix}.fasta', quality_annotations, module_number,
                    index_rep2=second_module_number, zip_it=zip_it, cutoff_after=cutoff_alignments)
    write_alignment(f'{motif_dir}/alignment_filtered_{suffix}.fasta', filt_primer, module_number,
                    index_rep2=second_module_number, zip_it=zip_it, cutoff_after=cutoff_alignments)
    write_alignment(f'{motif_dir}/alignment_filtered_left_{suffix}.fasta', [a for a in filt_primer if a.module_bases[0] > 0], module_number,
                    index_rep2=second_module_number, zip_it=zip_it, cutoff_after=cutoff_alignments)
    write_alignment(f'{motif_dir}/alignment_filtered_right_{suffix}.fasta', [a for a in filt_primer if a.module_bases[-1] > 0], module_number,
                    index_rep2=second_module_number, zip_it=zip_it, cutoff_after=cutoff_alignments, right_align=True)

    # write histogram image
    if second_module_number is not None:
        write_histogram_image2d(f'{motif_dir}/repetitions_{suffix}', quality_annotations + filt_primer, module_number, second_module_number,
                                motif_class.module_str(module_number), motif_class.module_str(second_module_number))
    else:
        write_histogram_image(f'{motif_dir}/repetitions_{suffix}', quality_annotations, filt_primer, module_number)
        write_profile(f'{motif_dir}/profile_{suffix}.txt', quality_annotations, index_rep=module_number)

    # write histogram txt files
    write_histogram(f'{motif_dir}/repetitions_{suffix}.txt', quality_annotations)
    write_histogram_nomenclature(f'{motif_dir}/nomenclatures_{suffix}.txt', quality_annotations, index_rep=module_number,
                                 index_rep2=second_module_number)
    write_histogram(f'{motif_dir}/repetitions_grey_{suffix}.txt', filt_primer)
    write_histogram_nomenclature(f'{motif_dir}/nomenclatures_grey_{suffix}.txt', filt_primer, index_rep=module_number,
                                 index_rep2=second_module_number)


def read_all_call(allcall_file: str) -> tuple[float, int | str, int | str, float, float, float, float, float, float] | None:
    """
    Read AllCall output and returns allele predictions and confidences.
    :param allcall_file: str - filename of AllCall output
    :return: tuple(9) - overall_confidence, allele numbers, and confidences for them, 0 for allele number if BG is the best
    """
    if not os.path.exists(allcall_file):
        return None

    with open(allcall_file) as f:
        lines = f.readlines()

    overall_conf = float(lines[0].strip().split()[-1].split('%')[0]) / 100

    def get_allele(line: str) -> tuple[int | str, float | str]:
        """
        Get allele number and its confidence from a line.
        :param line: str - a single line of AllCall output
        :return: (int/str, float) - allele number and its confidence
        """
        split = line.strip().split()
        num: int | str = split[0]
        conf: float | str = split[-1].split('%')[0].split(')')[0]
        try:
            conf = float(conf) / 100
        except ValueError:
            pass
        try:
            num = int(num)
        except ValueError:
            pass
        return num, conf

    def get_probability(line: str) -> float:
        perc = line.strip().split()[-1][:-1]
        return float(perc)

    a1, c1 = get_allele(lines[1])
    a2, c2 = get_allele(lines[2])
    c3 = get_probability(lines[3])
    c4 = get_probability(lines[4])
    c5 = get_probability(lines[5])
    c6 = get_probability(lines[6])

    return overall_conf, a1, a2, c1, c2, c3, c4, c5, c6


def custom_format(template: str, **kwargs) -> str:
    """
    Custom format of strings for only those that we provide
    :param template: str - string to format
    :param kwargs: dict - dictionary of strings to replace
    :return: str - formatted string
    """
    for k, v in kwargs.items():
        template = template.replace('{%s}' % k, v)

    return template


def find_file(filename: str, include_gzip: bool = False) -> typing.Optional[str]:
    """
    Find if we have a file with the provided filename.
    :param filename: str - file name
    :param include_gzip: bool - try the gzipped suffix
    :return: str/None - filename if exists, None if not
    """
    if os.path.exists(filename):
        return filename
    gzipped = filename + '.gz'
    if include_gzip and os.path.exists(gzipped):
        return gzipped
    return None


def generate_nomenclatures(filename: str, motif: Motif | None = None, nomenclature_limit: int = -1) -> list[str]:
    """
    Generate nomenclature string lines from nomenclature file. Maximally generate nomenclature_limit lines.
    :param filename: str - file name of the nomenclature file
    :param motif: Motif - motif class
    :param nomenclature_limit: int - limit how many nomenclatures to generate
    :return: list[str] - array of nomenclature lines
    """
    if not os.path.exists(filename):
        return []

    with open(filename, 'r') as noms:
        lines = []
        for line in noms:
            if line == '' or line is None:
                break

            line_split = line.split('\t')
            motif_parts = [f'<td>{s}</td>' for s in line_split[1:]]
            ref = f'{motif.chrom}:g.{motif.start}_{motif.end}' if motif is not None else ''
            nom_row = report.html_templates.nomenclature_string.format(count=line_split[0] + 'x', ref=ref, parts='\n    '.join(motif_parts))
            lines.append(nom_row)

            # end?
            if 0 < nomenclature_limit <= len(lines):
                break

    return lines


def write_report(motifs: list[Motif], result_table: pd.DataFrame, post_filter: PostFilter, report_dir: str, nomenclature_limit: int = 5) -> None:
    """
    Generate and write a report.
    :param motifs: list[Motif] - list of motifs to write results
    :param result_table: pd.DataFrame - table of results
    :param post_filter: PostFilter - post-filter arguments
    :param report_dir: str - dir name for reports
    :param nomenclature_limit: int - number of lines from nomenclature.txt to print
    :return: None
    """
    # tsv file with table:
    columns_to_rename = {'motif_name': 'Motif', 'motif_nomenclature': 'Nomenclature', 'motif_sequence': 'Sequence',
                         'repetition_index': 'Repetition index', 'confidence': 'Overall confidence', 'allele1': 'Allele 1 prediction',
                         'conf_allele1': 'Allele 1 confidence', 'reads_a1': 'Allele 1 reads', 'allele2': 'Allele 2 prediction',
                         'conf_allele2': 'Allele 2 confidence', 'reads_a2': 'Allele 2 reads', 'quality_reads': 'Reads (full)',
                         'one_primer_reads': 'Reads (partial)', 'conf_background_all': 'Both Background prob.',
                         'conf_background': 'One Background prob.', 'conf_extended_all': 'Background Expanded prob.',
                         'conf_extended': 'One Expanded prob.', 'indels': 'Indels per read', 'mismatches': 'Mismatches per read',
                         'indels_a1': 'Allele 1 indels', 'mismatches_a1': 'Allele 1 mismatches', 'indels_a2': 'Allele 2 indels',
                         'mismatches_a2': 'Allele 2 mismatches'}
    result_table = result_table[columns_to_rename.keys()]

    # merge all_profiles:
    all_profiles = f'{report_dir}/all_profiles.txt'
    all_true = f'{report_dir}/all_profiles.true'

    mcs = {}
    ms = {}
    rows = {}
    alignments = {}
    mcs_static = []
    ms_static = []
    rows_static = []

    with open(all_profiles, 'w') as pf, open(all_true, 'w') as tf:
        for motif in motifs:
            seq = motif.modules_str(include_flanks=True)
            for _, result in result_table[result_table['motif_name'] == motif.name].iterrows():

                # adjust helper variables
                phasing = '_' in str(result['repetition_index'])
                suffix = result['repetition_index']

                # read files
                rep_file = find_file(f'{report_dir}/{motif.dir_name()}/repetitions_{suffix}.json')
                pcol_file = find_file(f'{report_dir}/{motif.dir_name()}/pcolor_{suffix}.json')
                align_file = find_file(f'{report_dir}/{motif.dir_name()}/alignment_{suffix}.fasta', include_gzip=True)
                filt_align_file = find_file(f'{report_dir}/{motif.dir_name()}/alignment_filtered_{suffix}.fasta', include_gzip=True)
                filt_left_file = find_file(f'{report_dir}/{motif.dir_name()}/alignment_filtered_left_{suffix}.fasta', include_gzip=True)
                filt_right_file = find_file(f'{report_dir}/{motif.dir_name()}/alignment_filtered_right_{suffix}.fasta', include_gzip=True)
                nomenclature_lines = generate_nomenclatures(f'{report_dir}/{motif.dir_name()}/nomenclatures_{suffix}.txt', motif, nomenclature_limit)

                # generate rows of tables and images
                row = report.html_templates.generate_row(seq, result, post_filter)
                rows[motif.name] = rows.get(motif.name, []) + [row]
                rows_static.append(row)

                # add the tables
                mc, m, a = report.html_templates.generate_motifb64(seq, result, rep_file, pcol_file, align_file, filt_align_file, filt_left_file,
                                                                   filt_right_file, nomenclature_lines, post_filter)
                if motif.name in mcs:
                    ms[motif.name].append(m)
                    alignments[motif.dir_name()][1].append(a[1])
                else:
                    mcs[motif.name] = mc
                    ms[motif.name] = [m]
                    alignments[motif.dir_name()] = (a[0], [a[1]])

                # and static tables
                mc, m, a = report.html_templates.generate_motifb64(seq, result, rep_file, pcol_file, align_file, filt_align_file, filt_left_file,
                                                                   filt_right_file, nomenclature_lines, post_filter, static=True)
                if mc not in mcs_static:
                    mcs_static.append(mc)
                ms_static.append(m)

                if not phasing:
                    # add to profiles
                    with open(f'{report_dir}/{motif.dir_name()}/profile_{suffix}.txt') as po:
                        line = po.readline()
                        pf.write(f'{motif.name}_{suffix}\t{line}\n')

                    # add to true
                    tf.write(f'{motif.name}_{suffix}\t{result["allele1"]}\t{result["allele2"]}\n')

    # load the report template
    script_dir = os.path.dirname(os.path.abspath(__file__))
    template = open(f'{script_dir}/report.html', 'r').read()
    template_static = open(f'{script_dir}/report.html', 'r').read()

    # fill sample name and version of software
    sample = os.path.basename(report_dir)
    version = open(f'{os.path.dirname(os.path.dirname(script_dir))}/version.txt', 'r').read().strip()
    template = custom_format(template, sample=sample, version=version)
    template_static = custom_format(template_static, sample=sample, version=version)
    tabs = []

    # save the report file
    with open(f'{report_dir}/report.html', 'w') as f:
        contents_table = report.html_templates.contents.format(table='\n'.join(sorted(mcs.values())))
        template = custom_format(template, motifs_content=contents_table + '\n' + report.html_templates.make_datatable_string)

        # generate content for each motif
        for motif in sorted(motifs, key=lambda x: x.name):
            motif_clean = re.sub(r'[^\w_]', '', motif.name.replace('/', '_'))
            nomenclature_lines = generate_nomenclatures(f'{report_dir}/{motif.dir_name()}/nomenclature.txt', motif, nomenclature_limit)
            tabs.append(report.html_templates.motif_summary.format(motif_id=motif_clean,
                                                                   nomenclatures='\n'.join(nomenclature_lines), table='\n'.join(rows[motif.name]),
                                                                   motifs='\n'.join(ms[motif.name])))

        f.write(custom_format(template, table='', motifs='\n'.join(tabs)))

    # save the static report file
    with open(f'{report_dir}/report_static.html', 'w') as f:
        contents_table = report.html_templates.contents.format(table='\n'.join(mcs_static))
        table = report.html_templates.motif_summary_static.format(table='\n'.join(rows_static))

        f.write(custom_format(template_static, motifs_content=contents_table,
                              table=table, motifs='\n'.join(ms_static)))

    # write alignments as html files
    for motif_dir_name in alignments.keys():
        template_alignments = open(f'{script_dir}/alignments.html', 'r').read()
        template_alignments = custom_format(template_alignments, sample=motif_dir_name, motif_desc=alignments[motif_dir_name][0])

        with open(f'{report_dir}/{motif_dir_name}/alignments.html', 'w') as f:
            f.write(custom_format(template_alignments, alignments='\n'.join(alignments[motif_dir_name][1])))

    # copy javascript libraries
    shutil.copy2(f'{script_dir}/msa.min.gz.js', f'{report_dir}/msa.min.gz.js')
    shutil.copy2(f'{script_dir}/plotly-2.14.0.min.js', f'{report_dir}/plotly-2.14.0.min.js')
    shutil.copy2(f'{script_dir}/jquery-3.6.1.min.js', f'{report_dir}/jquery-3.6.1.min.js')
    shutil.copy2(f'{script_dir}/datatables.min.js', f'{report_dir}/datatables.min.js')

    # save the table(s) - we already generate a table in main file TODO unify
    result_table.rename(columns=columns_to_rename).to_csv(f'{report_dir}/table.tsv', sep='\t')
