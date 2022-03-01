# GUNC and BUSCO snakemake filtering workflow

This is a snakemake workflow made to integrate GUNC and BUSCO tools and run them on large sets of MAGs!
## Dependencies
It needs installed and working snakemake. Follow this [link](https://snakemake.readthedocs.io/en/stable/getting_started/installation.html) for details.

## Config file
It  is important to keep the formatting of the config file (pay attention to two spaces at the beggining of every row!).

This is the meaning of the following parameters:
1) genomes - path to the directory containing MAGs
2) input_text_file - text file containing paths to each genome being analysed. This enables having complex directory structures
3) temporary_dir - important when running on the HPC. Usually local scratch storage that has high I/O speeds
4) directory_faa - directory where to store .faa and .fna output of the Prodigal tool
5) format1, format2 - enables having MAGs with different extensions (zipped or unzipped, .fasta/.fna...)
6) database_folder - path to the dir where to download both GUNC and BUSCO databases
7) gunc_db_type - choose between gdtb and progenomes
8) busco_parameters - set either '-l *lineage*' (lineage = official lineage from BUSCO docs), --auto-lineage-prok or cluster analysis (cluster analysis is a specific function you most likely don't need)
8) batch_size - how many genomes in one batch
9) all_genomes, busco_output, cluster_info - specific fields necessary only if busco_parameters == 'cluster_analysis'
## How to run it

1) Create a working dir
2) Run the workflow with appropriate arguments, e.g. on HPC:
```bash
snakemake --default-resources mem_mb=12000 time=300 -j 10 -s ../workflow/Snakefile --profile cluster
```
Consult [snakemake docs](https://snakemake.readthedocs.io/en/stable/executing/cli.html) for more details about command-line interface.

## Contributing
This workflow was primarily made for internal use, thus not bug-proof. Please contact us if you need any specific help!
