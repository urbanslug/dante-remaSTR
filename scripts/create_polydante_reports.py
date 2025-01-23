from jinja2 import Environment, FileSystemLoader
import argparse
import json
import os
import pandas as pd  # change to polars?
import sys
import shutil
from typing import cast
from scipy.stats import ks_2samp  # type: ignore

# import pprint
# pp = pprint.PrettyPrinter()
# pp.pprint(data)


def load_arguments() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        formatter_class=argparse.RawDescriptionHelpFormatter,
        description="Generate reports comparing two groups as defined in polydante"
    )

    arg = parser.add_argument
    arg("-a", "--cases", nargs="+", help='Path to data.json files (possibly many) of samples with disease.')
    arg("-b", "--controls", nargs="+", help='Path to data.json files (possibly many) of healthy samples.')
    arg("-m", "--snv-map", help="Path to SNV map")
    arg("-o", "--output-dir", default="./results",
        help='Path to directory where the output will be stored. Default=<input_dir>/result_files')

    return parser.parse_args()


def append_data(filenames: list[str], group: str, motif2seq: dict[str, str]) -> list[tuple]:
    result = []
    for filename in filenames:
        print(f"Loading {filename}")
        with open(filename, "r") as f:
            data = json.load(f)
            for motif in data["motifs"]:
                str_id, snps = motif["motif_id"].split("-")
                snps2 = snps.split("_")
                assert len(motif["modules"]) == 1
                motif2seq[str_id] = motif["modules"][0]["sequence"]

                for snp in snps2:
                    data_blob = (
                        motif["modules"][0]["allele_1"],
                        motif["modules"][0]["allele_2"],
                        motif["modules"][0]["stats"],
                        motif["modules"][0]["reads_spanning"],
                        motif["modules"][0]["reads_flanking"],
                        motif["motif_stats"]
                    )
                    result.append((snp, str_id, data["sample"], group, data_blob))
    return result


def load_snv_map(filename: str) -> dict[str, str]:
    with open(filename) as f:
        lines = f.readlines()

    snv_id2rs_id = {}
    for line in lines[1:]:
        tmp = line.split()
        snv_id2rs_id[tmp[0]] = tmp[1]
    return snv_id2rs_id


def generate_table(df_str: pd.DataFrame) -> list[dict]:
    rows = []
    for idx, row in df_str.iterrows():
        _, _, sample, group, data = row
        tmp = {
            "sample": sample,
            "group": group,
            "a1_pred": data[0][0],
            "a1_conf": data[0][1],
            "a1_indel": data[0][2],
            "a1_mismatch": data[0][3],
            "a1_read": data[0][4],
            "a2_pred": data[1][0],
            "a2_conf": data[1][1],
            "a2_indel": data[1][2],
            "a2_mismatch": data[1][3],
            "a2_read": data[1][4],
            "conf": data[2][0],
            "indel": data[2][1],
            "mismatch": data[2][2],
            "read_full": data[3],
            "read_part": data[4]
        }
        rows.append(tmp)
    return rows


def allele_num(x: int | str) -> int:
    if x == "E":
        return -3
    if x == "X":
        return -2
    if x == "B":
        return -1
    return x  # type: ignore


def generate_heatmap_z(df_str: pd.DataFrame, max_allele: int):
    z1: list[list[int]] = [[0] * (max_allele + 1 + 2) for _ in range(max_allele + 1 + 2)]
    for _, row in df_str.iterrows():
        a1 = allele_num(row["data"][0][0])
        a2 = allele_num(row["data"][1][0])

        if a1 >= 0 and a2 >= 0:
            z1[a1 + 1][a2 + 1] += 1
        else:
            if a1 == -1 and a2 == -1:
                z1[0][0] += 1
                continue
            if a1 == -1:
                z1[0][a2 + 1] += 1
                continue
            if a1 == -3 and a2 == -3:
                z1[max_allele + 2][max_allele + 2] += 1
                continue
            if a2 == -3:
                z1[a1 + 1][max_allele + 2] += 1
                continue
            raise ValueError(f"{(a1, a2)} has unexpected value.")

    z2 = cast(list[list[int | None]], z1)
    for i in range(len(z1)):
        for j in range(i):
            assert z2[i][j] == 0, f"z[{i}][{j}] is non-zero"
            z2[i][j] = None
    return z2


def get_min_allele(df_str: pd.DataFrame) -> int:
    a1_preds = [allele_num(row["data"][0][0]) for _, row in df_str.iterrows()]
    a2_preds = [allele_num(row["data"][1][0]) for _, row in df_str.iterrows()]
    alleles = [x for x in a1_preds + a2_preds if x >= 0]
    return min(alleles)


def get_max_allele(df_str: pd.DataFrame) -> int:
    a1_preds = [allele_num(row["data"][0][0]) for _, row in df_str.iterrows()]
    a2_preds = [allele_num(row["data"][1][0]) for _, row in df_str.iterrows()]
    MIN_GS = 8  # minimal graph size, histograms/heatmaps with less columns look weird
    return max(a1_preds + a2_preds + [MIN_GS])
    # return max(a1_preds + a2_preds)  # TODO?


def generate_heatmap(df_str: pd.DataFrame) -> dict:
    max_allele = get_max_allele(df_str)
    # TODO?
    # min_allele = get_min_allele(df_str)
    # print(min_allele, max_allele)

    z_case = generate_heatmap_z(df_str.loc[df_str["group"] == "case"], max_allele)
    z_control = generate_heatmap_z(df_str.loc[df_str["group"] == "control"], max_allele)
    tickvals = list(range(max_allele + 3))
    ticktext = ["B"] + list(range(max_allele + 1)) + ["E"]
    xlim = max_allele + 1.5

    result = {
        "case": {
            "z": z_case,
            "tickvals": tickvals,
            "ticktext": ticktext,
            "xlim": xlim
        },
        "control": {
            "z": z_control,
            "tickvals": tickvals,
            "ticktext": ticktext,
            "xlim": xlim
        }
    }
    return result


def collect_histogram(df_str: pd.DataFrame, a_idx: int, max_allele: int) -> list[int]:
    result = [0] * (max_allele + 3)
    for _, row in df_str.iterrows():
        allele: int | str = row["data"][a_idx][0]
        if allele == 'E':
            result[-1] += 1
            continue
        if allele == 'B':
            result[0] += 1
            continue
        if allele == 'X':
            raise NotImplementedError

        allele = int(allele)
        result[allele + 1] += 1
    return result


def generate_histogram(df_str: pd.DataFrame) -> dict:
    max_allele = get_max_allele(df_str)

    df_case = df_str[df_str["group"] == "case"]
    zip_case = zip(collect_histogram(df_case, 0, max_allele), collect_histogram(df_case, 1, max_allele))
    y_case = [a + b for a, b in zip_case]

    df_control = df_str[df_str["group"] == "control"]
    zip_control = zip(collect_histogram(df_control, 0, max_allele), collect_histogram(df_control, 1, max_allele))
    y_control = [a + b for a, b in zip_control]

    tickvals = list(range(max_allele + 3))
    ticktext = ["B"] + list(range(max_allele + 1)) + ["E"]

    result = {
        "case": {
            "y": y_case,
            "tickvals": tickvals,
            "ticktext": ticktext,
        },
        "control": {
            "y": y_control,
            "tickvals": tickvals,
            "ticktext": ticktext,
        }
    }
    return result


def compare_statistically(histogram_data: dict) -> dict:
    cases = []
    for i, n in enumerate(histogram_data["case"]["y"][1:]):
        cases += [i] * n

    controls = []
    for i, n in enumerate(histogram_data["control"]["y"][1:]):
        controls += [i] * n

    if len(cases) == 0 or len(controls) == 0:
        statistic, p_value = float('nan'), float('nan')
    else:
        statistic, p_value = ks_2samp(cases, controls)

    return {"statistic": statistic, "p_value": p_value}


def generate_data(snv_id: str, df_snv: pd.DataFrame, motif2seq: dict[str, str]) -> dict:
    str_list = []
    for (str_id,), df_str in df_snv.groupby(["str"]):
        str_id = str_id
        seq = motif2seq[str_id]
        rows = generate_table(df_str)
        heatmap_data = generate_heatmap(df_str)
        histogram_data = generate_histogram(df_str)
        stat_results = compare_statistically(histogram_data)
        if stat_results["p_value"] < 0.05:
            print(f"{snv_id} {str_id}: The two distributions are significantly different.")
        str_location = df_str["data"].iloc[0][5]

        str_data = {
            "str_id": str_id,
            "str_seq": seq,
            "str_location": str_location,
            "stat_results": stat_results,
            "rows": rows,
            "heatmap_data": heatmap_data,
            "histogram_data": histogram_data
        }
        str_list.append(str_data)

    data = {
        "snv_id": snv_id,
        "str_list": str_list
    }

    return data


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


def main(args: argparse.Namespace) -> None:
    os.makedirs(args.output_dir, exist_ok=True)

    snv_id2rs_id: dict[str, str] = load_snv_map(args.snv_map)
    motif2seq: dict[str, str] = {}
    cases = append_data(args.cases, "case", motif2seq)
    controls = append_data(args.controls, "control", motif2seq)
    df = pd.DataFrame.from_records(cases + controls, columns=["snv", "str", "sample", "group", "data"])

    for (snv_id,), df_snv in df.groupby(["snv"]):
        snv_id = snv_id2rs_id[snv_id]
        print(f"Generating {args.output_dir}/{snv_id}.html")
        data = generate_data(snv_id, df_snv, motif2seq)

        script_dir = os.path.dirname(sys.argv[0]) + "/../templates"
        env = Environment(loader=FileSystemLoader([script_dir]), trim_blocks=True, lstrip_blocks=True)
        template = env.get_template("polydante_report.html")
        output = template.render(data=data)
        with open(f"{args.output_dir}/{snv_id}.html", "w") as f:
            f.write(output)

    copy_includes(args.output_dir)


if __name__ == '__main__':
    args = load_arguments()
    main(args)
