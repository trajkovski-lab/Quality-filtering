import glob
import os
import sys
with open(snakemake.log[0], 'w') as f:
    sys.stdout = f
    for i in snakemake.input.working_dir:
        command_os = f"mkdir {snakemake.output}"
        os.system(command_os)
        glob_argument = f"{i}/*.gff"
        for path in glob.glob(glob_argument):
            sample_name = "".join(path.split("/")[-1].split(".")[:-1])
            command = f"any2fasta {path} > {snakemake.output}/{sample_name}.fasta"
            os.system(command)
        unzipped_format = snakemake.config['format1'].strip(".gz")
        glob_argument_fasta = f"{i}/*{unzipped_format}"
        if glob.glob(glob_argument_fasta):
            for path_fasta in glob.glob(glob_argument_fasta):
                sample_name = "".join(path_fasta.split("/")[-1][:-(len(unzipped_format))])+"fasta"
                new_path = "/".join(path_fasta.split("/")[:-1])+"/"+sample_name
                command_rename = f"mv {path_fasta} {new_path}"
                os.system(command_rename)
                command = f"mv {new_path} {snakemake.output}"
                os.system(command)
