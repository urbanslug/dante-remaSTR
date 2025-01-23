#!/usr/bin/env python3
from common import HGVSRecord, chr_to_num

input = "/dev/stdin"
with open(input) as f:
    lines = f.readlines()

output = []
for line in lines:
    disease_id, hgvs, desc = line.split("\t")
    hgvs_rec = HGVSRecord(hgvs)

    chrom = hgvs_rec.chrom
    start = hgvs_rec.start - 1
    end = hgvs_rec.end
    name = disease_id + "__" + str(hgvs_rec)
    new_line = f"{chrom}\t{start}\t{end}\t{name}"
    output.append((chr_to_num[chrom], start, new_line))

output.sort()
for item in output:
    print(item[2])
