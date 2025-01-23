#!/usr/bin/env python3

from Bio import SeqIO
import sys


class ExpHunterRecord:
    # CHROM POS ID REF ALT QUAL FILTER INFO FORMAT *sample_IDs
    def __init__(self, vcf_line: str):
        items = vcf_line.strip().split()
        self.chrom = items[0]
        self.pos = int(items[1])
        self.id = items[2]
        self.ref = items[3]
        self.alt = items[4]
        self.qual = items[5]
        self.filter = items[6]
        self.info: dict[str, str | bool] = process_info(items[7])
        self.format = items[8]
        # self.sample_IDs = list(map(lambda x: process_samples(items[8], x), items[9:]))
        self.sample_IDs = items[9:]

    def __repr__(self) -> str:
        info = info_str(self.info)
        samples = "\t".join(self.sample_IDs)
        # CHROM POS ID REF ALT QUAL FILTER INFO FORMAT *sample_IDs
        return "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}".format(
            self.chrom, self.pos, self.id, self.ref, self.alt,
            self.qual, self.filter, info, self.format, samples
        )


def normalize_ref_alt(ref: str, alt: str) -> tuple[str, str]:
    suffix = 0
    it = zip(reversed(ref), reversed(alt))
    for _ in range(min(len(ref), len(alt)) - 1):
        (a, b) = next(it)
        if a == b:
            suffix += 1

    return (ref[:len(ref) - suffix], alt[:len(alt) - suffix])


def process_info(item: str) -> dict[str, str | bool]:
    tmp: dict[str, str | bool] = dict(map(lambda x: info_to_pair(x), item.split(";")))
    return tmp


def info_to_pair(item: str) -> tuple[str, str | bool]:
    split = item.split("=")
    match len(split):
        case 1:
            return (split[0], True)
        case 2:
            return (split[0], split[1])
        case _:
            raise ValueError(f"INFO contains {item}")


def info_str(info: dict[str, str | bool]) -> str:
    info_str = []
    for key, val in info.items():
        if type(val) is bool:
            info_str.append(key)
        else:
            info_str.append(f"{key}={val}")
    return ";".join(info_str)


# vcf = "./EH_HG002_results_norm.vcf"
# ref_seq = "../data/grch38_decoy.fa"
vcf = "/dev/stdin"
ref_seq = sys.argv[1]

with open(vcf) as f:
    lines = f.readlines()

with open(ref_seq) as f:
    tmp = list(SeqIO.parse(f, "fasta"))
    fa_records = dict([(x.id, x) for x in tmp])

# print header
header_lines = list(filter(lambda x: x.startswith("##"), lines))
print(header_lines[0], end="")
for line in header_lines[1:]:
    if line.startswith("##ALT"):
        continue
    if line.startswith("##contig"):  # ExpansionHunter forgets contig length...
        continue
    print(line, end="")
print('##INFO=<ID=SVLEN,Number=1,Type=Integer,Description="Alternate length - Reference length">')
print('##contig=<ID=chr1,length=248956422>')
print('##contig=<ID=chr2,length=242193529>')
print('##contig=<ID=chr3,length=198295559>')
print('##contig=<ID=chr4,length=190214555>')
print('##contig=<ID=chr5,length=181538259>')
print('##contig=<ID=chr6,length=170805979>')
print('##contig=<ID=chr7,length=159345973>')
print('##contig=<ID=chr8,length=145138636>')
print('##contig=<ID=chr9,length=138394717>')
print('##contig=<ID=chr10,length=133797422>')
print('##contig=<ID=chr11,length=135086622>')
print('##contig=<ID=chr12,length=133275309>')
print('##contig=<ID=chr13,length=114364328>')
print('##contig=<ID=chr14,length=107043718>')
print('##contig=<ID=chr15,length=101991189>')
print('##contig=<ID=chr16,length=90338345>')
print('##contig=<ID=chr17,length=83257441>')
print('##contig=<ID=chr18,length=80373285>')
print('##contig=<ID=chr19,length=58617616>')
print('##contig=<ID=chr20,length=64444167>')
print('##contig=<ID=chr21,length=46709983>')
print('##contig=<ID=chr22,length=50818468>')
print('##contig=<ID=chrX,length=156040895>')
print('##contig=<ID=chrY,length=57227415>')

# print records
record_lines = lines[len(header_lines):]
print(record_lines[0], end="")
for line in record_lines[1:]:
    record = ExpHunterRecord(line)
    if record.alt == ".":
        continue

    # we need to take one letter to the left, because
    # otherwise we couldn't represent alt_copy_number == 0
    sequence = fa_records[record.chrom].seq
    prev_letter = sequence[record.pos - 1 - 1]  # 1-based to 0-based and left

    ref_copy_number = int(record.info["REF"])
    ref_seq = prev_letter + record.info["RU"] * ref_copy_number

    alt_copy_number = int(record.alt[len("<STR"):-1])
    alt_seq = prev_letter + record.info["RU"] * alt_copy_number

    svlen = len(alt_seq) - len(ref_seq)
    assert svlen != 0, "Oopsie daisy!"
    svtype = "INS" if svlen > 0 else "DEL"
    (ref_seq, alt_seq) = normalize_ref_alt(ref_seq, alt_seq)

    record.pos -= 1  # took one letter left
    record.ref = ref_seq
    record.alt = alt_seq
    record.info['SVLEN'] = str(svlen)
    record.info['SVTYPE'] = svtype
    del record.info['END']
    print(record)
