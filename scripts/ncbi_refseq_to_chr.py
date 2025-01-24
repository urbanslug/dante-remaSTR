#!/usr/bin/env python3
from common import HGVSRecord, chr_to_num

input = "/dev/stdin"
with open(input) as f:
    lines = f.readlines()

output = []
for line in lines[1:]:
    disease_id, hgvs = line.split("\t")
    hgvs_rec = HGVSRecord(hgvs)
    new_line = f"{disease_id}\t{hgvs_rec}\n"
    output.append((chr_to_num[hgvs_rec.chrom], hgvs_rec.start, new_line))

output.sort()
for item in output:
    print(item[2], end="")
