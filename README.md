# GUNC and BUSCO snakemake filtering workflow

This is a snakemake workflow made to integrate GUNC and BUSCO tools and run them on large sets of MAGs in batches!


## Dependencies
It needs installed and working snakemake. Follow this [link](https://snakemake.readthedocs.io/en/stable/getting_started/installation.html) for details.

## Config file

See also in [`config/default_config.yaml`](config/default_config.yaml)
```yaml
## It  is important to keep the formatting of the config file (pay attention to two spaces at the beggining of every row!).

## Input

# option A: text file  with absolute paths to genomes are being analyzed one per line 
input_text_file: genome_paths.txt 


# option B: all genomes in a folder (set input_text_file to "")
genomes: /data/genomes/ #path to genome folder


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




## How to run it

1) Create a working dir
2) Run the workflow with appropriate arguments, e.g. on HPC:
```bash
snakemake --default-resources mem_mb=12000 time=300 -j 10 -s ../workflow/Snakefile --profile cluster
```
Consult [snakemake docs](https://snakemake.readthedocs.io/en/stable/executing/cli.html) for more details about command-line interface.

## Cluster workflow

The cluster workflow is intended for **re-evaluation of genome quality estimation** on a per species basis. The commonly used tools for genome quality estimation, e.g. CheckM and BUSCO have the option to automatically choose the best reference set of marker genes for the given genome. This is very usefull, but can also induce some bias. For example, let's take a set of genomes expected to be from the same species (e.g. ANI <= 95%). Genomes in this set might be evaluated on strikingly different reference set because they were automatically selected by BUSCO. In the extreme case a very "bad" genome could have a better quality score than a bery "good" genome just because it was evaluated on the root marker gene set. 

To avoid this we came up with the idea of this workflow where we run the Workflow once on all species (cluster) representetives and then use the same busco lineage for all genomes in this cluster. This guarantees that all genomes in a cluster are compared with the same measure. It also accelerates the BUSCO run.

We made this workflow to re-filter the genomes from HumGut but the aproach can be generalized to any genomes that are clustered. 

TODO: For the cluster workflow you need 







## Contributing
This workflow was primarily made for internal use, thus not bug-proof. Please contact us if you need any specific help!
