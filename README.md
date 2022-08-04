# GUNC and BUSCO snakemake filtering workflow

This is a snakemake workflow made to integrate GUNC and BUSCO tools and run them on large sets of MAGs in batches!


## Dependencies
It needs installed and working snakemake. Follow this [link](https://snakemake.readthedocs.io/en/stable/getting_started/installation.html) for details.

## Config file

See also in [`config/default_config.yaml`](config/default_config.yaml)
```yaml
## It  is important to keep the formatting of the config file (pay attention to two spaces at the beggining of every row!).

## Input

# option A: text file  with absolute paths to genomes, one per line 
input_text_file: directory_structure.txt


# option B: all genomes in a folder (set input_text_file to "")
genomes: genome #path to a dir with genomes


# input format to be expected: enables having MAGs with different extensions (zipped or unzipped, .fasta/.fna, gff, anything that goes trough any2fasta)
format1: .fna.gz
format2: .fasta.gz

batch_size: 50 # how many genomes in one batch

temporary_dir: /tmp # important when running on the HPC. Usually local scratch storage that has high I/O speeds
directory_faa: faa_files # (sub)directory where to store .faa and .fna output of the Prodigal tool
database_folder: ../../resources/database/ # path to the dir where to download both GUNC and BUSCO databases

gunc_db_type: gtdb # choose between gtdb and progenomes

busco_parameters: --auto-lineage-prok #  set either '-l *lineage*' (lineage = official lineage from BUSCO docs), --auto-lineage-prok or cluster analysis (see below)
save_busco_output: False # If True, write the whole busco output in main output dir. Else, only tables containing quality metrics will be kept.

#Resources
threads: 8 # number of threads

mem_mb_prodigal: 12000 # memory in MB
time_min_prodigal: 30 # time in minutes

mem_mb_gunc: 24000 # memory in MB
time_min_gunc: 60 # time in minutes

mem_mb_busco: 12000 # memory in MB
time_min_busco: 180 # time in minutes


# specific fields necessary only if busco_parameters == 'cluster_analysis'
all_genomes: ../../All_genomes.tsv # table where cluster is specified for every genome in genomes/ dir
busco_output: ../../BUSCO_hq_reps_output.tsv #output for BUSCO representative run
cluster_info: cluster95 #Which cluster was used to create previous busco output

```




## How to test it
1) Clone this repo and make a test_run directory inside it:
```bash
git clone https://github.com/trajkovski-lab/Quality-filtering.git
cd Quality-filtering
mkdir test_run
```
2) Download and unzip test data from [this Zenodo link:](https://zenodo.org/record/6322941)
```bash
cd test_run
wget https://zenodo.org/record/6322941/files/genome.tar.gz
wget https://zenodo.org/record/6322941/files/directory_structure.txt
tar -xvzf genome.tar.gz
```
3) Use dry-run feature of snakemake to test the workflow and if everything is working, run it:
```bash
snakemake -s ../workflow/Snakefile -j 5 -c 4 -np
snakemake -s ../workflow/Snakefile -j 5 -c 4
```
Consult [snakemake docs](https://snakemake.readthedocs.io/en/stable/executing/cli.html) for more details about command-line interface.

## Cluster workflow

The cluster workflow is intended for **re-evaluation of genome quality estimation** on a per species basis. The commonly used tools for genome quality estimation, e.g. CheckM and BUSCO have the option to automatically choose the best reference set of marker genes for the given genome. This is very useful, but can also induce some bias. For example, let's take a set of genomes expected to be from the same species (e.g. ANI <= 95%). Genomes in this set might be evaluated on strikingly different reference set because they were automatically selected by BUSCO. In the extreme case a very "bad" genome could have a better quality score than a very "good" genome just because it was evaluated on the root marker gene set. 

To avoid this we came up with the idea of this workflow where we run the Workflow once on all species (cluster) representatives and then use the same BUSCO lineage for all genomes in this cluster. This guarantees that all genomes in a cluster are compared with the same measure. It also accelerates the BUSCO run.

We made this workflow to re-filter the genomes from HumGut but the approach can be generalised to any genomes that are clustered. 

TODO: For the cluster workflow you need 







## Contributing
This workflow was primarily made for internal use, thus not bug-proof. Please contact us if you need any specific help!
