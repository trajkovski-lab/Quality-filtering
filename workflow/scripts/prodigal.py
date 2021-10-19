from glob import glob
from snakemake.shell import shell
import os, sys
import logging, traceback

logging.basicConfig(
    filename= snakemake.log[0],
    level=logging.INFO,
    format="%(asctime)s %(message)s",
    datefmt="%Y-%m-%d %H:%M:%S",
)

logging.captureWarnings(True)


def handle_exception(exc_type, exc_value, exc_traceback):
    if issubclass(exc_type, KeyboardInterrupt):
        sys.__excepthook__(exc_type, exc_value, exc_traceback)
        return

    logging.error(
        "".join(
            [
                "Uncaught exception: ",
                *traceback.format_exception(exc_type, exc_value, exc_traceback),
            ]
        )
    )

# Install exception handler
sys.excepthook = handle_exception


from multiprocessing import Pool

def run_prodigal(input_fasta):
    # run prodigal on one sample
    # stderr get raised and saved to log file
    sample_name = os.path.splitext(os.path.basename(input_fasta))[0]
    logging.info(f"Run prodigal on sample {sample_name}")
    shell(f"prodigal {snakemake.params.parameters} -i {input_fasta} -a {snakemake.output}/{sample_name}.faa > /dev/null")




def run_multiple_prodigal(input_dir, out_dir, log, threads,extension='.fasta'):

    # read input list
    genomes_fastas = glob(os.path.join(input_dir, "*"+extension))
    assert len(genomes_fastas)>0
    os.makedirs(out_dir, exist_ok=True)

    pool = Pool(int(threads))
    pool.map( run_prodigal,genomes_fastas)


for input in snakemake.input:
    lineage = "".join("".join("".join(input.split("/")[-1]).split(".")[0]).split("_")[4:])
    with open(input,'r') as path_file:
        os.system(f"mkdir -p /scratch/{snakemake.config['intermediate_dir']}/results_any2fasta/{snakemake.wildcards.counter}-{snakemake.wildcards.lineage}")
        paths = path_file.readlines()
        for path in paths:
            path_command = path.strip("\n")
            unzipped_sample = "".join(path_command.split("/")[-1]).rstrip(".gz")
            sample_name = ".".join(unzipped_sample.split(".")[:-1])
            if "gff" in unzipped_sample or ".gz" in path_command:
                shell(f"any2fasta {path_command} > /scratch/{snakemake.config['intermediate_dir']}/results_any2fasta/{snakemake.wildcards.counter}-{snakemake.wildcards.lineage}/{sample_name}.fasta")
            else:
                dst = f"/scratch/{snakemake.config['intermediate_dir']}/results_any2fasta/{snakemake.wildcards.counter}-{snakemake.wildcards.lineage}/{sample_name}.fasta"
                shell(f"ln -s {path_command} {dst}")
            os.system(f"mkdir -p {snakemake.output}")
            run_prodigal(f"/scratch/{snakemake.config['intermediate_dir']}/results_any2fasta/{snakemake.wildcards.counter}-{snakemake.wildcards.lineage}/{sample_name}.fasta")


'''
os.system(f"mkdir -p {snakemake.output}")
run_multiple_prodigal(
    f"{snakemake.config['intermediate_dir']}/results_any2fasta/{snakemake.wildcards.counter}-{snakemake.wildcards.lineage}/", snakemake.output[0], snakemake.log[0], int(snakemake.threads))
'''
