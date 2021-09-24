import os
import glob

for i in snakemake.input.fasta_dir:
    command_os = f"mkdir {snakemake.output}"
    os.system(command_os)
    glob_argument = f"{i}/*.fasta"
    for path in glob.glob(glob_argument):
        sample_name = path.split("/")[-1].split(".")[0]
        command = f"prodigal -q -p meta -i {path} -a {snakemake.output}/{sample_name}.faa &> {log}"
        os.system(command)
