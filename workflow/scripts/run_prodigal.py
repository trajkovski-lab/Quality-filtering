import glob
import os
import sys
from snakemake.shell import shell


os.makedirs(snakemake.output[0])

for path in glob.glob( f"{snakemake.input.fasta_dir}/*.fasta"):
    #I don't understand this rule
    sample_name = "".join(path.split("/")[-1].split(".")[:-1])
    shell(f"prodigal {snakemake.params.parameters} -i {path} -a {snakemake.output}/{sample_name}.faa &> {log}")
