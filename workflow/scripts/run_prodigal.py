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


import glob
import os
import sys
from snakemake.shell import shell


os.makedirs(snakemake.output[0])





for path in glob.glob( f"{snakemake.input.fasta_dir}/*.fasta"):
    sample_name = os.path.splitext(os.path.basename(path))[0]
    logging.info(f"Run prodigal on sample {sample_name}")

    shell(f"prodigal {snakemake.params.parameters} -i {path} -a {snakemake.output[0]}/{sample_name}.faa  ")
