import re
import sys


ncbi_refseq_to_chr: dict[str, str] = {
    "NC_000001": "chr1",  "NC_000002": "chr2",  "NC_000003": "chr3",
    "NC_000004": "chr4",  "NC_000005": "chr5",  "NC_000006": "chr6",
    "NC_000007": "chr7",  "NC_000008": "chr8",  "NC_000009": "chr9",
    "NC_000010": "chr10", "NC_000011": "chr11", "NC_000012": "chr12",
    "NC_000013": "chr13", "NC_000014": "chr14", "NC_000015": "chr15",
    "NC_000016": "chr16", "NC_000017": "chr17", "NC_000018": "chr18",
    "NC_000019": "chr19", "NC_000020": "chr20", "NC_000021": "chr21",
    "NC_000022": "chr22", "NC_000023": "chrX"
}


chr_to_num: dict[str, int] = {
    "chr1" :  1, "chr2" :  2, "chr3" :  3, "chr4" :  4, "chr5" :  5, "chr6" :  6,
    "chr7" :  7, "chr8" :  8, "chr9" :  9, "chr10": 10, "chr11": 11, "chr12": 12,
    "chr13": 13, "chr14": 14, "chr15": 15, "chr16": 16, "chr17": 17, "chr18": 18,
    "chr19": 19, "chr20": 20, "chr21": 21, "chr22": 22, "chrX" : 23, "chrY" : 24
}


class HGVSRecord:
    # NC_000006.12:g.45422680_45422802AGC[24]CGG[17]
    # chr1:g.45422680_45422802AGC[24]CGG[17]
    def __init__(self, nomenclature: str) -> None:
        m1 = re.match(r"(NC_[0-9]+\.[0-9]+|chr[0-9XY]+):g\.([0-9]+)_([0-9]+)(.+)", nomenclature)
        if m1 is None:
            raise ValueError(f"{nomenclature} does not have a correct format.")
        seqid, start, end, lstruct = m1.groups()

        # remap any id to chr
        if seqid.startswith("NC"):
            seqid = ncbi_refseq_to_chr[seqid.split(".")[0]]

        # parse motifs
        m2 = re.findall(r"([A-Z]+\[[0-9]+\])", lstruct)
        locus_struct: list[tuple[str, int]] = []
        for unit in m2:
            m3 = re.match(r"([A-Z]+)\[([0-9]+)\]", unit)
            if m3 is None:
                raise ValueError(f"Unit {unit} cannot be parsed.")
            seq, occ = m3.groups()
            locus_struct.append((str(seq), int(occ)))

        self.chrom = seqid
        self.start = int(start)
        self.end = int(end)
        self.units = locus_struct

    def __str__(self) -> str:
        units: list[str] = []
        for unit in self.units:
            units.append(f"{unit[0]}[{unit[1]}]")
        motifs = "".join(units)
        return f"{self.chrom}:g.{self.start}_{self.end}{motifs}"

    def get_loci_intervals(self) -> list[tuple[str, int, int]]:
        """
        chr1:g.100_119GGCGCGGAGC[2] -> [(chr1, 100, 119)]
        chr1:g.680_802AGC[24]CGG[17] -> [(chr1, 680, 751), (chr1, 752, 802)]
        """
        result = []
        start = self.start
        for seq, occ in self.units:
            end = start + len(seq) * occ
            result.append((self.chrom, start, end - 1))
            start = end
        if self.end == end:
            print(f"Warning: Motif {self} is inconsistent.", file=sys.stderr)
        return result
