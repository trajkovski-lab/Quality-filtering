import glob
import os
import sys
mkdir_command = f"mkdir {snakemake.output}"
os.system(mkdir_command)
with open(snakemake.input.lists, 'r') as path_file:
    paths = path_file.readlines()
    for path in paths:
        path_command = path.strip("\n")
        unzipped_sample = "".join(path_command.split("/")[-1]).rstrip(".gz")
        sample_name = ".".join(unzipped_sample.split(".")[:-1])
        if "gff" in unzipped_sample or ".gz" in path_command:
            command = f"any2fasta {path_command} > {snakemake.output}/{sample_name}.fasta"
            os.system(command)
        else:
            dst = f"{snakemake.output}/{sample_name}.fasta"
            command1 = f"ln -s {path_command} {dst}"
            os.system(command1)
