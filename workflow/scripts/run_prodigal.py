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

# Start of the scripts


from glob import glob
import os
import sys
from snakemake.shell import shell
from multiprocessing import Pool

def run_prodigal(input_fasta):
    # run prodigal on one sample
    # stderr get raised and saved to log file
    sample_name = os.path.splitext(os.path.basename(input_fasta))[0]
    logging.info(f"Run prodigal on sample {sample_name}")
    shell(f"prodigal {snakemake.params.parameters} -i {input_fasta} -a {snakemake.output[0]}/{sample_name}.faa > /dev/null")

    


def run_multiple_prodigal(input_dir, out_dir, log, threads,extension='.fasta'):

    # read input list
    genomes_fastas = glob(os.path.join(input_dir, "*"+extension))
    assert len(genomes_fastas)>0
    os.makedirs(out_dir, exist_ok=True)

    pool = Pool(threads)
    pool.map( run_prodigal,genomes_fastas)





if __name__ == "__main__":


    run_multiple_prodigal(
        snakemake.input.fasta_dir,
        snakemake.output[0],
        snakemake.log[0],
        int(snakemake.threads),
    )
