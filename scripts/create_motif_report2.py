from collections import defaultdict
from typing import TypeAlias
import argparse
import json
import textwrap


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


def collect_jsons(inputs: list[str]) -> dict[str, list[tuple[str, dict]]]:
    motif_dict: dict[str, list] = defaultdict(list)
    for filename in inputs:
        print(f"Reading {filename}")
        with open(filename, "r") as f:
            data = json.load(f)

        for motif in data["motifs"]:
            motif_dict[motif["motif_id"]].append((data["sample"], motif))

    return motif_dict


def main() -> None:
    args = load_arguments()
    motif_dict = collect_jsons(args.inputs)

    tmp_dict = {}
    for k, v in motif_dict.items():
        if k == "ALS" or k == "DM2":
            tmp_dict[k] = v
    print(tmp_dict)


if __name__ == '__main__':
    main()
