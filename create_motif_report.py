from collections import defaultdict
from jinja2 import Environment, FileSystemLoader
from typing import TypeAlias
import argparse
import json
import os
import textwrap
import numpy as np


HistReadCounts: TypeAlias = tuple[list[int], list[int], list[int]]
Row: TypeAlias = tuple[str, str | int, str, str | int, str, str, int, int, HistReadCounts]
HeatmapData: TypeAlias = tuple
HistData: TypeAlias = list[int]
PopulationData: TypeAlias = tuple[str, str, int, HeatmapData, HistData, list[Row]]

MIN_GS = 10  # minimal graph size, histograms/heatmaps with less columns look weird


def load_arguments() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        formatter_class=argparse.RawDescriptionHelpFormatter,
        description=textwrap.dedent("")
    )

    arg = parser.add_argument
    arg("-i", "--inputs", nargs="+", help='Path to data.json files (possibly many).')
    arg("-o", "--output_dir", default="./results",
        help='Path to directory where the output will be stored. Default=<input_dir>/result_files')
    arg("--report_every", default=5, type=int,
        help='Specify how often a progress message should be printed (default=5)')
    arg('-q', '--quiet', action='store_true', help='Don\'t print any progress messages')

    return parser.parse_args()


def collect_alleles(hist_data: list[int], bg_num: int, rows: list[Row], a_idx: int) -> tuple[list[int], int]:
    for row in rows:
        allele: int | str = row[a_idx]  # type: ignore  # if the index is correct, this is correct type
        if allele == 'E':
            hist_data[-1] += 1
            continue
        if allele == 'B':
            bg_num += 1
            continue
        if allele == 'X':
            raise NotImplementedError

        allele = int(allele)
        hist_data[allele] += 1
    return hist_data, bg_num


def collect_heatmap_data(max_allele: int, pairs: list[tuple[int | str, int | str]]) -> HeatmapData:
    z: list[list[int]] = [[0] * (max_allele + 1 + 2) for _ in range(max_allele + 1 + 2)]

    for pair in pairs:
        a1 = allele_num(pair[0])
        a2 = allele_num(pair[1])

        if a1 >= 0 and a2 >= 0:
            z[a1 + 1][a2 + 1] += 1
        else:
            if a1 == -1 and a2 == -1:
                z[0][0] += 1
                continue
            if a1 == -1:
                z[0][a2 + 1] += 1
                continue
            if a1 == -3 and a2 == -3:
                z[max_allele + 2][max_allele + 2] += 1
                continue
            if a2 == -3:
                z[a1 + 1][max_allele + 2] += 1
                continue
            raise ValueError(f"{pair} has unexpected value.")

    z_new: list[list[int | None]] = z  # type: ignore
    for i in range(len(z)):
        for j in range(i):
            assert z_new[i][j] == 0, f"z[{i}][{j}] is non-zero"
            z_new[i][j] = None
    tickvals = list(range(max_allele + 3))
    ticktext = ["B"] + list(range(max_allele + 1)) + ["E"]
    xlim = max_allele + 1.5

    heatmap_data = (z_new, tickvals, ticktext, xlim)
    return heatmap_data


def main(args: argparse.Namespace) -> None:
    loci_dict: dict[str, list[Row]] = defaultdict(list[Row])
    loci_to_seq: dict = {}
    for filename in args.inputs:
        print(f"Reading {filename}")
        with open(filename, "r") as f:
            data = json.load(f)

        (sample, _, _, motifs) = data
        for motif in motifs:
            (_, _, loci, _) = motif
            for locus in loci:
                (locus_id, sequence, _, allele1, allele2, stats, spanning, flanking, graph_data) = locus
                (a1_pred, a1_conf, _, _, _) = allele1
                (a2_pred, a2_conf, _, _, _) = allele2
                (confidence, _, _) = stats
                (histogram_data, _, _) = graph_data

                row = (sample, a1_pred, a1_conf, a2_pred, a2_conf, confidence, spanning, flanking, histogram_data)
                loci_dict[locus_id].append(row)

                if locus_id not in loci_to_seq:
                    loci_to_seq[locus_id] = sequence
                else:
                    assert loci_to_seq[locus_id] == sequence, f"{loci_to_seq[locus_id]} != {sequence}: Mixed motif DB?"

    print("Collecting done. Starting to compute aggregates.")
    loci_data: dict[str, PopulationData] = {}
    for locus_id, rows in loci_dict.items():
        sequence = loci_to_seq[locus_id]
        allele_pairs = [(a1, a2) for _, a1, _, a2, _, _, _, _, _ in rows]

        # ----------------------------------------------------------------------
        a1_max = max(allele_num(row[1]) for row in rows)
        a2_max = max(allele_num(row[3]) for row in rows)
        max_allele = max(a1_max, a2_max, MIN_GS)
        # print(locus_id, max_allele)

        # collect bg_num and hist_data
        bg_num = 0
        # hist_data = [0] * (max_allele + 2)  # +1 for zero indexing, +1 for expansion
        hist_data = [0] * (40 + 2)  # +1 for zero indexing, +1 for expansion
        hist_data, bg_num = collect_alleles(hist_data, bg_num, rows, 3)
        hist_data, bg_num = collect_alleles(hist_data, bg_num, rows, 1)

        heatmap_data = collect_heatmap_data(max_allele, allele_pairs)
        # ----------------------------------------------------------------------

        loci_data[locus_id] = (locus_id, sequence, bg_num, heatmap_data, hist_data, rows)

    print("Aggregates done. Creating HTMLs.")
    script_dir = os.path.dirname(__file__)
    os.makedirs(f"{args.output_dir}/aggregate", exist_ok=True)
    for locus_id, data in loci_data.items():
        env = Environment(loader=FileSystemLoader([script_dir]), trim_blocks=True, lstrip_blocks=True)
        template = env.get_template("./dante_remastr_standalone_templates/population_report.html")
        output = template.render(main_data=data)
        with open(f"{args.output_dir}/aggregate/{locus_id}.html", "w") as f:
            f.write(output)


def allele_num(x: int | str) -> int:
    if x == "E":
        return -3
    if x == "X":
        return -2
    if x == "B":
        return -1
    return x  # type: ignore


if __name__ == '__main__':
    args = load_arguments()
    main(args)
