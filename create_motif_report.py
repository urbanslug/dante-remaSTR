from itertools import chain
import argparse
import os
import re
import json
import textwrap
import numpy as np
from jinja2 import Environment, FileSystemLoader
from typing import TextIO, Iterator, TypeAlias, Any
from collections import defaultdict


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


def main(args: argparse.Namespace) -> None:
    loci_dict: dict[str, list[Row]] = defaultdict(list[Row])
    loci_to_seq: dict = {}
    for filename in args.inputs:
        print(filename)
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
        # collect bg_num
        bg_num = 0
        # collect heatmap_data
        heatmap_data = (1, )
        # collect hist_data
        hist_data = (1, )

        loci_data[locus_id] = (locus_id, sequence, bg_num, heatmap_data, hist_data, rows)

    script_dir = os.path.dirname(__file__)
    os.mkdir(f"{args.output_dir}")
    os.mkdir(f"{args.output_dir}/aggregate")
    for locus_id, data in loci_data.items():
        env = Environment(loader=FileSystemLoader([script_dir]), trim_blocks=True, lstrip_blocks=True)
        template = env.get_template("./dante_remastr_standalone_templates/population_report.html")
        output = template.render(main_data=data)
        with open(f"{args.output_dir}/aggregate/{locus_id}.html", "w") as f:
            f.write(output)


HistReadCounts: TypeAlias = tuple[list[int], list[int], list[int]]
Row: TypeAlias = tuple[str, str | int, str, str | int, str, str, int, int, HistReadCounts]
HeatmapData: TypeAlias = tuple
HistData: TypeAlias = tuple
PopulationData: TypeAlias = tuple[str, str, int, HeatmapData, HistData, list[Row]]
if __name__ == '__main__':
    args = load_arguments()
    main(args)
