# GUNC and BUSCO snakemake filtering workflow

This is a snakemake workflow made to integrate GUNC and BUSCO tools and run them on large sets of MAGs!
## Dependencies
It needs installed and working snakemake. Follow this [link](https://snakemake.readthedocs.io/en/stable/getting_started/installation.html) for details.

## Config file

See also in [`config/default_config.yaml`](config/default_config.yaml)
```yaml
## It  is important to keep the formatting of the config file (pay attention to two spaces at the beggining of every row!).

## Input
genomes: /data/genomes/ #path to the base genome folder
input_text_file: rikenellaceae.txt # txt file with paths to genomes are being analyzed

# input format to be expected: enables having MAGs with different extensions (zipped or unzipped, .fasta/.fna, gff, anything that goes trough any2fasta)
format1: .fasta.gz
format2: .fna.gz

batch_size: 50 # how many genomes in one batch

temporary_dir: /tmp # important when running on the HPC. Usually local scratch storage that has high I/O speeds
directory_faa: faa_files # (sub)directory where to store .faa and .fna output of the Prodigal tool
database_folder: ../../resources/database/ # path to the dir where to download both GUNC and BUSCO databases

gunc_db_type: gtdb # choose between gdtb and progenomes

busco_parameters: --auto-lineage-prok #  set either '-l *lineage*' (lineage = official lineage from BUSCO docs), --auto-lineage-prok or cluster analysis (see below)


# specific fields necessary only if busco_parameters == 'cluster_analysis'
all_genomes: ../../All_genomes.tsv #table where cluster is specified for every genome in genomes/ dir
busco_output: ../../BUSCO_hq_reps_output.tsv #output for BUSCO representative run
cluster_info: cluster95 #Which cluster was used to create previous busco output
```

Example of a `input_text_file`:




## How to run it

1) Create a working dir
2) Run the workflow with appropriate arguments, e.g. on HPC:
```bash
snakemake --default-resources mem_mb=12000 time=300 -j 10 -s ../workflow/Snakefile --profile cluster
```
Consult [snakemake docs](https://snakemake.readthedocs.io/en/stable/executing/cli.html) for more details about command-line interface.

## Contributing
This workflow was primarily made for internal use, thus not bug-proof. Please contact us if you need any specific help!
