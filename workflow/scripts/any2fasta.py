import glob
import os
import sys
from snakemake.shell import shell


os.makedirs(snakemake.output[0])

with open(snakemake.input.lists, 'r') as path_file:
    paths = path_file.readlines()

for path in paths:
    path_command = path.strip("\n")
    unzipped_sample = "".join(path_command.split("/")[-1]).rstrip(".gz")
    sample_name = ".".join(unzipped_sample.split(".")[:-1])

    if "gff" in unzipped_sample or ".gz" in path_command:
        #run any2fasta
        shell(f"any2fasta {path_command} > {snakemake.output}/{sample_name}.fasta")

    else:
        # simply create smlink
        os.symlink(os.path.abspath(path_command), f"{snakemake.output}/{sample_name}.fasta")
