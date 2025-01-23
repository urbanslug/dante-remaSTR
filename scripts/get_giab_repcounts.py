#!/usr/bin/env python3
"""Calculate ground truth from giab sample."""

import gzip
from common import HGVSRecord

REFERENCE = "./../inputs/GRCh38_decoy.fa"
GIAB_VCF = "./../cache/HG002_GRCh38_diseases_nosnp_v1.0.1.vcf.gz"
MOTIFS = "./../cache/diseases_predominant_dante.tsv"
SLOP = 30


def get_svlen(chrom: str, start: int, end: int, lines: list[str]) -> tuple[int, int]:
    a1_len = end + 1 - start  # end + 1 because [start, end] is closed interval
    a2_len = end + 1 - start

    for line in lines:
        # CHROM POS ID REF ALT QUAL FILTER INFO FORMAT *sample_IDs
        vcf_chrom, p, _, ref, alt, _, _, _, _format, *samples = line.split()
        pos = int(p)
        gt = samples[0].split(":")[0].split("|")

        if vcf_chrom != chrom or pos < start or pos > end:
            continue

        # print(vcf_chrom, pos, ref, alt, samples)

        diff = len(alt) - len(ref)
        if gt[0] == "1":
            a1_len += diff
        if gt[1] == "1":
            a2_len += diff

    return a1_len, a2_len


def main() -> None:
    with gzip.open(GIAB_VCF, "rb") as f1:
        tmp = map(lambda x: x.decode("utf-8"), f1.readlines())
        giab_lines = list(filter(lambda x: not x.startswith("#"), tmp))

    with open(MOTIFS, "r", encoding="utf-8") as f2:
        motif_lines = f2.readlines()

    for x in motif_lines:
        name, nom = x.split()
        nomenclature = HGVSRecord(nom)
        loci_intervals = nomenclature.get_loci_intervals()
        n = len(loci_intervals)
        for i, (chrom, start, end) in enumerate(loci_intervals):
            # print(f"{name}\t{nomenclature}\t{chrom}\t{start}\t{end}")
            unit_len: int = len(nomenclature.units[i][0])
            if i == 0 and i == n - 1:  # for locus which is starting and ending at the same time
                a1_len, a2_len = get_svlen(chrom, start - SLOP, end + SLOP, giab_lines)
                a1 = (a1_len - 2 * SLOP) / unit_len
                a2 = (a2_len - 2 * SLOP) / unit_len
            elif i == 0:
                a1_len, a2_len = get_svlen(chrom, start - SLOP, end, giab_lines)
                a1 = (a1_len - SLOP) / unit_len
                a2 = (a2_len - SLOP) / unit_len
            elif i == n - 1:
                a1_len, a2_len = get_svlen(chrom, start, end + SLOP, giab_lines)
                a1 = (a1_len - SLOP) / unit_len
                a2 = (a2_len - SLOP) / unit_len
            else:
                a1_len, a2_len = get_svlen(chrom, start, end, giab_lines)
                a1 = (a1_len) / unit_len
                a2 = (a2_len) / unit_len
            # print(f"{a1}\t{a2}")
            smaller = min(a1, a2)
            bigger = max(a1, a2)
            print(f"{name}\t{i}\t{nomenclature}\t{a1:.2f}\t{a2:.2f}\t{smaller:.2f}\t{bigger:.2f}")
        # print()


if __name__ == "__main__":
    main()
