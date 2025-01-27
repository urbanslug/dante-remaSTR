from collections import defaultdict, Counter
from typing import TypeAlias
from jinja2 import Environment, FileSystemLoader
from typing import Any
import argparse
import json
import textwrap
import os
import sys
import shutil
import pandas as pd


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


def copy_includes(output_dir: str) -> None:
    include_dir = os.path.dirname(sys.argv[0]) + "/../includes"
    os.makedirs(f'{output_dir}/includes', exist_ok=True)
    shutil.copy2(f'{include_dir}/msa.min.gz.js',            f'{output_dir}/includes/msa.min.gz.js')
    shutil.copy2(f'{include_dir}/plotly-2.14.0.min.js',     f'{output_dir}/includes/plotly-2.14.0.min.js')
    shutil.copy2(f'{include_dir}/jquery-3.6.1.min.js',      f'{output_dir}/includes/jquery-3.6.1.min.js')
    shutil.copy2(f'{include_dir}/datatables.min.js',        f'{output_dir}/includes/datatables.min.js')
    shutil.copy2(f'{include_dir}/styles.css',               f'{output_dir}/includes/styles.css')
    shutil.copy2(f'{include_dir}/w3.css',                   f'{output_dir}/includes/w3.css')
    shutil.copy2(f'{include_dir}/jquery.dataTables.css',    f'{output_dir}/includes/jquery.dataTables.css')


def convert_modules(modules: list) -> list[str]:
    return [f"{seq}[{n}]" for seq, n in modules]


def create_main_histrogram(df: pd.DataFrame) -> dict:
    counter = Counter(list(df["nomenclature1"]) + list(df["nomenclature2"]))

    data = [(v, k, i) for i, (k, v) in enumerate(counter.items())]
    data.sort()

    result = {
        "tickvals": [x[2] for x in data],
        "x": [x[0] for x in data],
        "ticktext": [x[1] for x in data]
    }
    return result


def main() -> None:
    args = load_arguments()
    motif_dict = collect_jsons(args.inputs)

    print("Aggregates done. Creating HTMLs.")
    for motif, v in motif_dict.items():
        if motif != "ALS" and motif != "DM2":
            continue

        print(f"Creating report fo motif {motif}")
        # print(json.dumps(v, indent=4))
        data: dict[str, Any] = {}

        data["sequence"] = convert_modules(v[0][1]["motif_stats"]["modules"])
        # create main table
        data_tmp1 = []
        for sample, tmp_data in v:
            data_tmp1.append((
                sample,
                tmp_data["phased_seqs"]["nomenclature1"],
                tmp_data["phased_seqs"]["errors1"],
                tmp_data["phased_seqs"]["nomenclature2"],
                tmp_data["phased_seqs"]["errors2"]
            ))

        columns = ["sample", "nomenclature1", "warnings1", "nomenclature2", "warnings2"]
        df = pd.DataFrame.from_records(data_tmp1, columns=columns)
        data["main_table"] = json.loads(df.to_json(orient="records"))
        # print(df.to_json(orient="records", indent=4))

        data["main_histogram"] = create_main_histrogram(df)
        # print(json.dumps(tmp_dict, indent=4))

        template_dir = os.path.dirname(__file__) + "/../templates"
        os.makedirs(f"{args.output_dir}", exist_ok=True)
        env = Environment(loader=FileSystemLoader([template_dir]), trim_blocks=True, lstrip_blocks=True)
        template = env.get_template("population_report2.html")
        output = template.render(data=data)
        print(f"Writting {args.output_dir}/{motif}.html")
        with open(f"{args.output_dir}/{motif}.html", "w") as f:
            f.write(output)

    copy_includes(args.output_dir)


if __name__ == '__main__':
    main()
