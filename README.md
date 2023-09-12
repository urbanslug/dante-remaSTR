# Dante (remastered)
Dante ("Da Amazing NucleoTide Exposer") is an algorithm for genotyping STR alleles based on NGS reads originating from the STR locus of interest. The method accounts for natural deviations from the expected sequence, such as variation in the repeat count, sequencing errors, ambiguous bases, and complex loci containing several motifs.

Reported figures provide an evidence for an expanded allele, which is too long to be captured by a single NGS read as well as for allelic single point mutations, small insertions, and deletions that may be relevant for a diagnostic evaluation.


## Getting Started

These instructions will get you a copy of the project up and running on your local machine. 

### Installation

Dante is developed and tested on Python 3.11 with libraries determined in the `conda_env.yaml`. Create the conda environment as following (although 'conda' works, it is slow, and we recommend to use 'mamba'):
```bash
mamba env create -n dante -f conda_env.yaml
```

If you want to use `create_motif_report.py` script to aggregate reports by motif, you will need to install additional `beautifulsoup4` and `natsort` dependencies:

```bash
mamba install beautifulsoup4 natsort
```

## Creating motif reports

To collect results from multiple reports and group them by motif, we provide a script `create_motif_report.py`. It has one required and one optional argument: `input dir` and `output dir`. Input dir is a path to folder with reports from samples. There can be other files as well, the script filters out all files that don't match `*.html`. Output dir is a path to directory, where the motif reports will be generated. If the path isn't specified, the script generates reports to `example/motif_report`.

After running Dante on example dataset, you may run this script as:

```bash
python create_motif_report.py example/report [example/motif_report]
```

Script generates separate report for each unique motif in Dante reports. Script uses the file name of the report to differentiate the source of result in table, so it is recommended to rename the reports before starting this script, so no two reports have the same name.

