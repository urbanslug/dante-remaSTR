#!/usr/bin/env python3
import sys

from common import HGVSRecord


def get_ref_region(locus_struct, locus_position):
    chrom, start, end = locus_position
    pos = start
    result = []
    for seq, occ in locus_struct:
        motif_len = len(seq) * occ
        # chrX:p-p+2 represents 3 bases => subtract 1 from the end
        # result.append(f'"{chrom}:{pos}-{pos + motif_len - 1}"')
        # but ExpansionHunter is shifted also at the start
        result.append(f'"{chrom}:{pos - 1}-{pos + motif_len - 1}"')
        pos += motif_len
    if pos - 1 != end:
        print(f"Warning: inconsistent {chrom}:{start}-{end} and {start=} {locus_struct=}", file=sys.stderr)
    return ', '.join(result)


def get_json_obj(data: tuple[str, list[tuple[str, int]], tuple[str, int, int]]) -> str:  # noqa: E501
    locus_id, locus_struct, locus_position = data
    locus_structure = ''.join([f"({seq})*" for seq, _ in locus_struct])

    reference_region = get_ref_region(locus_struct, locus_position)
    variant_id = ', '.join([f'"{locus_id}_{seq}"' for seq, _ in locus_struct])
    variant_type = ', '.join(['"Repeat"'] * len(locus_struct))
    if len(locus_struct) > 1:
        reference_region = f"[{reference_region}]"
        variant_id = f"[{variant_id}]"
        variant_type = f"[{variant_type}]"

    return f"""
    {{
        "LocusId": "{locus_id}",
        "LocusStructure": "{locus_structure}",
        "ReferenceRegion": {reference_region},
        "VariantId": {variant_id},
        "VariantType": {variant_type}
    }}"""


def parse_line(line: str) -> tuple[str, list[tuple[str, int]], tuple[str, int, int]]:  # noqa: E501
    locus_id, hgvs, _ = line.split("\t")
    hgvs_rec = HGVSRecord(hgvs)
    locus_position = (hgvs_rec.chrom, hgvs_rec.start, hgvs_rec.end)
    return locus_id, hgvs_rec.units, locus_position


infile = "/dev/stdin"
with open(infile) as f:
    lines = f.readlines()

    data1 = map(parse_line, lines)
    data2 = map(get_json_obj, data1)
    data3 = ",".join(data2)
    print(f"[{data3}\n]")
