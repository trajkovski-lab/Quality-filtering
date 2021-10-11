import glob
import os

for i in snakemake.input.fasta_dir:
    command1 = f"mkdir {snakemake.output}"
    os.system(command1)
    glob_argument = f"{i}/*.faa"
    for path in glob.glob(glob_argument):
        sample_name = "".join(path.split("/")[-1].split(".")[:-1])
        command = f"busco -i {path} --auto-lineage-prok -m protein -o {sample_name} --out_path {snakemake.output} -c {snakemake.threads} &> {snakemake.log}"
        os.system(command)
