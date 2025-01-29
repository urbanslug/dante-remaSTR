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


def create_ticks(df: pd.DataFrame) -> list[str]:
    nom1 = set(df.loc[:, ["nomenclature1_len", "nomenclature1_join"]]
                 .itertuples(index=False, name=None))

    nom2 = set(df.loc[:, ["nomenclature2_len", "nomenclature2_join"]]
                 .itertuples(index=False, name=None))

    data = list(nom1.union(nom2))
    data.sort()
    ticktext = [x[1] for x in data]
    return ticktext


def create_main_histrogram(df: pd.DataFrame, ticks: list[str]) -> dict:
    data = [0] * len(ticks)
    for _, row in df.iterrows():
        idx = ticks.index(row["nomenclature1_join"])
        data[idx] += 1
        idx = ticks.index(row["nomenclature2_join"])
        data[idx] += 1

    result = {
        "tickvals": list(range(len(ticks))),
        "ticktext": ticks,
        "x": data
    }
    return result


def create_main_heatmap(df: pd.DataFrame, ticks: list[str]) -> dict:
    data: list[list[int | None]] = [[0 for _ in range(len(ticks))] for _ in range(len(ticks))]
    for _, row in df.iterrows():
        idx1 = ticks.index(row["nomenclature1_join"])
        idx2 = ticks.index(row["nomenclature2_join"])
        idx1, idx2 = min(idx1, idx2), max(idx1, idx2)
        data[idx1][idx2] += 1  # type: ignore

    for i in range(len(ticks)):
        for j in range(i + 1, len(ticks)):
            data[j][i] = None

    result = {
        "tickvals": list(range(len(ticks))),
        "ticktext": ticks,
        "z": data
    }
    return result


def generate_df(v):
    data_tmp1 = []
    for sample, tmp_data in v:
        data_tmp1.append((
            sample,
            tuple(tmp_data["phased_seqs"]["nomenclature1"]),
            tmp_data["phased_seqs"]["nomenclature1_len"],
            tmp_data["phased_seqs"]["errors1"],
            tuple(tmp_data["phased_seqs"]["nomenclature2"]),
            tmp_data["phased_seqs"]["nomenclature2_len"],
            tmp_data["phased_seqs"]["errors2"]
        ))

    columns = ["sample", "nomenclature1", "nomenclature1_len", "warnings1", "nomenclature2", "nomenclature2_len", "warnings2"]
    df = pd.DataFrame.from_records(data_tmp1, columns=columns)
    df["nomenclature1_join"] = df["nomenclature1"].apply(lambda x: "".join(x))
    df["nomenclature2_join"] = df["nomenclature2"].apply(lambda x: "".join(x))
    return df


def main() -> None:
    args = load_arguments()
    motif_dict = collect_jsons(args.inputs)

    print("Aggregates done. Creating HTMLs.")
    for motif, v in motif_dict.items():
        if motif != "ALS" and motif != "DM2":
            continue

        print(f"Creating report fo motif {motif}")

        sequence = convert_modules(v[0][1]["motif_stats"]["modules"])
        df = generate_df(v)
        ticks = create_ticks(df)

        data: dict[str, Any] = {}
        data["sequence"] = sequence
        data["main_table"] = json.loads(df.to_json(orient="records"))
        data["main_histogram"] = create_main_histrogram(df, ticks)
        data["main_heatmap"] = create_main_heatmap(df, ticks)
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
