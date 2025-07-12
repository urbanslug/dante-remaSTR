# Dante (remastered)
Dante ("Da Amazing NucleoTide Exposer") is an algorithm for genotyping STR alleles based on NGS reads originating 
from the STR locus of interest. The method accounts for natural deviations from the expected sequence, such as variation 
in the repeat count, sequencing errors, ambiguous bases, and complex loci containing several motifs.

Reported figures provide an evidence for an expanded allele, which is too long to be captured by a single NGS read 
as well as for allelic single point mutations, small insertions, and deletions that may be relevant for a diagnostic 
evaluation.

## Project structure

The whole algorithm is divided into two parts: 
1. annotation of mapped reads - `remastr` tool in Rust language due to parallelism and speed
    - managed in a separate repository
    - input: BAM file with mapped reads, txt file with HGVS nomenclature, fasta file with reference
    - output: tsv file with annotated reads, columns are: 
        - motif: motif nomenclature
        - read_sn: number of read for each motif
        - read_id: read identifier as in BAM
        - mate_order: information about read pairing: 1 - left read, 2 - right read, 0 - unpaired
        - read: read sequence
        - reference: aligned annotated motif in the read
        - modules: aligned numbers of motive modules 
        - log_likelihood: log likelihood of the annotation
2. genotyping of annotated reads, visualization and reporting - `dante_remastr` tool in Python language
    - input: tsv file with annotated reads (output of `remastr` tool)
    - output: tsv file with genotyped read motives and confidences:
        - motif_name: motif nomenclature
        - motif_sequence: motif STR (part of the motif)
        - chromosome: chromosome of the motif
        - start: start location of the motif
        - end: end location of the motif 
        - allele1: number of repeats of the STR (or B - background, E - expanded) on first allele
        - allele2: number of repeats of the STR (or B - background, E - expanded) on second allele 
        - confidence: confidence of the prediction in percents (counted as proportion of the likelihood of the predicted state versus all)
        - conf_allele1: confidence of the prediction of the first allele
        - conf_allele2: confidence of the prediction of the second allele
        - quality_reads: number of quality reads for this STR with both primers
        - one_primer_reads: number of quality reads for this STR with exactly one primer
        - filtered_reads: filtered reads mapped on this location (either no primer, modules are not consecutive, too many errors, or not 
                    enough of repetitions or bases in modules)
        - conf_background: confidence of the background prediction for both alleles
        - conf_background_all: confidence of the background prediction for at least one allele
        - conf_extended: confidence of the expanded allele prediction for both alleles
        - conf_extended_all: confidence of the expanded allele prediction for at least one allele

## Getting Started

These instructions will get you a copy of the project up and running on your local machine. 

### Acquisition of repository

This repository can be cloned with:

```bash
git clone https://github.com/marcelTBI/dante-remaSTR.git
```

### Conda environment

Dante is developed with libraries determined in the `conda_env.yaml`. Create the conda 
environment as following (although `conda` works, it is slow, and we recommend to use `mamba`):

```bash
mamba env create -n env_dante -f conda_env.yaml
```

### Download remastr subtool
```bash
git clone https://gitlab.com/andrejbalaz/remastr.git
```

#### Compilation and running

The `remastr` Rust part needs to be compiled. Navigate to the submodule and run compilation:

```bash
cd remastr/dante_cli
conda activate dante_remastr
cargo build --release  # rust needs to be installed (https://www.rust-lang.org/tools/install)
```

The compiled release version of the `dante_cli` tool should appear in the `remastr/target/release` directory after
successful compilation. 

The `dante_remastr` Python part does not need to be compiled.  

Run both parts as follows:

```bash
./remastr/target/release/dante_cli -b <reads_bam> -m <STRSet_file.tsv> -o remastr_output.tsv
python ./dante-remaSTR/dante_remastr_standalone.py -i remastr_output.tsv -o remastr_output -v
```

### Visualisation

The tool is able to provide comprehensive visualisations and user report in the form of a `.html` file of the results 
with the switch `--verbose` or `-v`. However, this is feasible only for up to about 100 motifs. You can run the smaller 
example file with additional visualisations as:

```bash
# go to the example directory
cd example
# activate the environment
conda activate dante_remastr
# run the genotyping and visualization
python dante_remastr_standalone.py -i remastr_output.tsv -o remastr_output -v
```
