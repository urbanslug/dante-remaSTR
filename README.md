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
    - managed in a separate repository and imported as a git submodule  
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
git clone --recurse-submodules git@github.com:marcelTBI/dante-remaSTR.git
```

The `--recurse-submodules` argument is needed because of the `remastr` submodule managed in a separate repositorry.

### Conda environment

Dante is developed with libraries determined in the `conda_env.yaml`. Create the conda 
environment as following (although `conda` works, it is slow, and we recommend to use `mamba`):

```bash
mamba env create -n dante_remastr -f conda_env.yaml
```

If you want to use `create_motif_report.py` script to aggregate reports by motif, you will need to install 
additional `beautifulsoup4` and `natsort` dependencies:

```bash
mamba install beautifulsoup4==4.12.2 natsort==7.1.1
```

The conda environment file defines environment for both parts of the project. 

#### Compilation and running

The `remastr` Rust part needs to be compiled. Navigate to the submodule and run compilation:

```bash
cd remastr
conda activate dante_remastr
cargo build --release
```

The compiled release version of the `remastr` tool should appear in the `remastr/target/release` directory after
successful compilation. 

The `dante_remastr` Python part does not need to be compiled.  

Run both parts as follows:

```bash
remastr/target/release/remastr -f <reference_fasta> -m <nomenclature_txt> -b <reads_bam> -o <remastr_output_tsv>
conda activate dante_remastr
python dante_remastr.py [optional_arguments] < <remastr_output_tsv> > <dante_output_tsv>
```

To list the optional arguments for the Python part, run:

```bash
python dante_remastr.py --help
```

Or simply use the provided pipeline script (try `python remastr_pipeline.py --help` for argument listing):

```bash
python remastr_pipeline.py <bams_files> <nomenclature_file> <output_dir>
```
 
### Example of the whole pipeline

We provide two examples (small and big) to try the program. 
To run the smaller one do:

```bash
# go to the example directory
cd example
# extract the chromosome X reference
gzip -d chromosomeX.fna.gz  
# run the annotation
../remastr/target/release/remastr -f chromosomeX.fna -n small_HGVS.txt -b small.bam -o small_generated.tsv
# activate the environment
conda activate dante_remastr
# run the genotyping
python ../dante_remastr.py --output-dir small_results < small_generated.tsv > small_res_generated.tsv
```

Note: the annotation will produce error messages in form `Read ... does not have pair information.`. This is expected. 

The results should be the same as those provided in the example directory, i.e. these commands should return no changes:
```bash
diff small_generated.tsv small.tsv
diff small_res_generated.tsv small_res.tsv
```

For the bigger example, we provide only data for the genotyping part (providing the BAM file and reference file would be too bulky for the repository). 
The bigger example is gzipped, which should be specified in the input arguments as the following: 

```bash
# go to the example directory
cd example
conda activate dante_remastr
# run the genotyping
python dante_remastr.py --output-dir big_results --input-gzipped < big.tsv.gz > big_res.tsv
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
python dante_remastr.py --output-dir small_results --verbose < small.tsv > small_res.tsv
```

## Creating motif reports

To collect results from multiple reports and group them by motif, we provide a script `create_motif_report.py`. 
It has one required and one optional argument: `input dir` and `output dir`. Input dir is a path to folder with 
reports from samples. There can be other files as well, the script filters out all files that don't match `*.html`.
Output dir is a path to directory, where the motif reports will be generated. If the path isn't specified, 
the script generates reports to `example/motif_report`.

After running Dante on example dataset, you may run this script as:

```bash
python create_motif_report.py example/report [example/motif_report]
```

Script generates separate report for each unique motif in Dante reports. Script uses the file name of the 
report to differentiate the source of result in table, so it is recommended to rename the reports before 
starting this script, so no two reports have the same name.

