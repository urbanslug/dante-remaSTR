# Dante (remastered)

<!-- markdown-toc start - Don't edit this section. Run M-x markdown-toc-refresh-toc -->
**Table of Contents**

- [Dante (remastered)](#dante-remastered)
- [Usage](#usage)
  - [Project structure](#project-structure)
    - [Output from `remastr` ](#output-from-remastr)
    - [Output from `dante_remastr` ](#output-from-dante_remastr)
  - [dante_cli](#dante_cli)
  - [The standalone dante python script](#the-standalone-dante-python-script)
- [Installation](#installation)
  - [Installation from Source](#installation-from-source)
    - [Clone the Repository](#clone-the-repository)
    - [Compile dante_cli](#compile-dante_cli)

<!-- markdown-toc end -->


Dante ("Da Amazing NucleoTide Exposer") is an algorithm for genotyping Short Tandem Repeat (STR) alleles based on Next Generation Sequencing (NGS) reads originating
from the STR locus of interest. The method accounts for natural deviations from the expected sequence, such as variation
in the repeat count, sequencing errors, ambiguous bases, and complex loci containing several motifs.

Reported figures provide a evidence for an expanded allele, which is too long to be captured by a single NGS read
as well as for allelic single point mutations, small insertions, and deletions that may be relevant for a diagnostic
evaluation.


# Usage

## Project structure

dante is divided into two parts: `remastr` and `dante_remastr`.

|     |     |     |
| ----| ----| ----|
| name| `remastr` | `dante_remastr` |
| task | annotation of mapped reads | genotyping of annotated reads, visualization and reporting |
| source | managed in a separate repository (https://gitlab.com/andrejbalaz/remastr) | from this repo  (https://github.com/marcelTBI/dante-remaSTR) |
| input |BAM file with mapped reads, txt file with HGVS nomenclature, fasta file with reference | tsv file with annotated reads (the input for `dante_remastr` is the output of `remastr` tool) |
|output | tsv file with annotated reads, columns are: | tsv file with genotyped read motives and confidences: |

### Output from `remastr` 

`remastr` tsv output has the columns as follows:
|column | purpose |
|--|-- |
| motif | motif nomenclature|
| read_sn | number of read for each motif|
| read_id | read identifier as in BAM|
|  mate_order | information about read pairing: 1 - left read, 2 - right read, 0 - unpaired|
|   read | read sequence |
|   reference | aligned annotated motif in the read |
|   modules | aligned numbers of motive modules |
|  log_likelihood | log likelihood of the annotation |

### Output from `dante_remastr` 

`dante_remastr` tsv output has the columns as follows:
|column | purpose |
|--|-- |
| motif_name | motif nomenclature|
| motif_sequence | motif STR (part of the motif)| 
| chromosome | chromosome of the motif |
| start | start location of the motif |
| end | end location of the motif |
| allele1 | number of repeats of the STR (or B - background, E - expanded) on first allele |
| allele2 | number of repeats of the STR (or B - background, E - expanded) on second allele |
| confidence | confidence of the prediction in percents (counted as proportion of the likelihood of the predicted state versus all)|
| conf_allele1 | confidence of the prediction of the first allele |
| conf_allele2 | confidence of the prediction of the second allele |
| quality_reads | number of quality reads for this STR with both primers
| one_primer_reads | number of quality reads for this STR with exactly one primer|
| filtered_reads | filtered reads mapped on this location (either no primer, modules are not consecutive, too many errors, or not enough of repetitions or bases in modules)|
| conf_background | confidence of the background prediction for both alleles|
| conf_background_all | confidence of the background prediction for at least one allele |
| conf_extended | confidence of the expanded allele prediction for both alleles|
| conf_extended_all | confidence of the expanded allele prediction for at least one allele |


## dante_cli

if `dante_cli` was both compiled and installed (e.g., using the command [`cargo install --path .`](#compile-dante_cli)) it should be accessible directly from your shell. If `dante_cli` was compiled but not installed, you can find the binary in the `remastr/target/release` directory.

`dante_cli` expects a sorted BAM file sorted by coordinates. 
If your BAM file is not sorted, you can use [samtools](https://github.com/samtools/samtools) to sort it by running the following command:
```bash
samtools sort <file_name>.bam -o <file_name>.sorted.bam
```

If you installed `dante_cli` locally, you can run it directly.
```bash
dante_cli -b <reads_bam>.sorted.bam -m <STRSet_file>.tsv -o remastr_output.tsv
```

If dante is not installed you have to specify the path to the binary when you run it as follows:
```bash
./remastr/target/release/dante_cli -b <reads_bam>.sorted.bam -m <STRSet_file>.tsv -o remastr_output.tsv
```

## The standalone dante python script

The project root `dante-remaSTR` contains a python script
[`dante_remastr_standalone.py`](dante_remastr_standalone.py) that we use to run the genotyping and visualise
the results from dante (described above in [Project structure](#project-structure)).
The script provides comprehensive visualisations and user report in the form of a `.html` file of the results
with the switch `--verbose` or `-v`. However, this is feasible only for up to about 100 motifs.

If the `env_dante` virtual environment, created earlier, is not active, activate it.
```bash
conda activate env_dante
```

Run the following to generate an analysis our results.
The output will be stored in the directory specified with `-o`.
```bash
python dante_remastr_standalone.py -i remastr_output.tsv -o remastr_output -v
```

You can run the smaller example file with additional visualisations using
```bash
python dante_remastr_standalone.py -i remastr_output.tsv -o remastr_output -v
```


# Installation

Currently, `dante` can only be installed from source.

## Installation from Source


Follow these instructions to get a copy of the project up and running on your local machine.


### Clone the Repository

Clone the `dante_remaSTR` repository and the `remastr` submodule
```bash
git clone --recursive https://github.com/marcelTBI/dante-remaSTR.git
```

Update the `remastr` submodule to the latest version (optional)
```bash
git submodule update --remote --merge
git commit -m "Update remastr submodule $(date +"%y%m%d%H%M%S")"
```

Dante depends on a set of libraries, which are managed through [mamba](https://anaconda.org/conda-forge/mamba), a C++ reimplementation of the conda package manager. These dependencies are specified in the [`conda_env.yaml`](./conda_env.yaml) file. 
If `mamba` isn't already installed on your system follow the installation instructions [here](https://mamba.readthedocs.io/en/latest/installation/mamba-installation.html) to install it.
Once mamba is installed, you can create the environment using the following command (while `conda` will also work, `mamba` is recommended for speed):

```bash
mamba env create -n env_dante -f conda_env.yaml
```

### Compile `dante_cli`

The `dante_remastr` python part does not need to be compiled.
However, part of `remastr` is written in rust and thus needs to be compiled.
If rust isn't already installed on your system follow the installation instructions [here](https://www.rust-lang.org/tools/install) to install it.


Navigate to the `dante_cli` subtool directory and compile the subtool as follows:
```bash
cd remastr/dante_cli
conda activate env_dante
```

Compile `dante_cli`

```bash
cargo build --release
```

To install `dante_cli` locally run
```bash
cargo install --path .
```

For more on cargo
```bash
cargo help # general cargo documentation
cargo help install # for documentation on installation
cargo help uninstall # for documentation on removing installed packages
```
