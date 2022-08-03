import glob
import os
import pandas as pd
from snakemake.shell import shell
import logging, traceback

logging.basicConfig(
    filename=snakemake.log[0],
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


def busco_lineage(path, lineage, specific_name, db_path, out_dir, threads, log, copy=False):
    shell(f"busco -i {path} -l {lineage} -f -q -m protein --offline -o {specific_name}--download_path {db_path} --out_path {out_dir} -c {threads} &> {log}")
    if copy:
        os.makedirs(snakemake.output.busco_output, exist_ok=True)
        shell(f"cp {out_dir}/{specific_name}/batch_summary.txt {snakemake.output.busco_output}")


def busco_autolineage(path, specific_name, db_path, out_dir, threads, log, copy=False):
    shell(f"busco -i {path} --auto-lineage-prok -m protein -f -q --offline -o {specific_name} --download_path {db_path} --out_path {out_dir} -c {threads} &> {log}")
    if copy:
        os.makedirs(snakemake.output.busco_output, exist_ok=True)
        shell(f"cp {out_dir}/{specific_name}/batch_summary.txt {snakemake.output.busco_output}")


temporary_dir = snakemake.config['temporary_dir']
path = snakemake.input.fasta_dir
db_path = snakemake.config["database_folder"]
specific_name = f"{snakemake.wildcards.counter}-{snakemake.wildcards.lineage}"

save_output = snakemake.config['save_busco_output']
if save_output:
    out_dir = "quality_filtering/intermediate_results/busco_output"
    if snakemake.config["busco_parameters"] == "cluster_analysis":
        lineage = "".join("".join(path.split("/")[-2]).split("-")[-1])
        busco_lineage(path, lineage, specific_name, db_path, out_dir, snakemake.threads, snakemake.log)

    elif "auto-lineage-prok" in snakemake.config["busco_parameters"]:
        busco_autolineage(path, specific_name, db_path, out_dir, snakemake.threads, snakemake.log)

    else:
        lineage = snakemake.config["busco_parameters"]
        busco_lineage(path, lineage, specific_name, db_path, out_dir, snakemake.threads, snakemake.log)


elif save_output == False:
    out_dir = f"{temporary_dir}/quality_filtering/intermediate_results/busco_output"
    if snakemake.config["busco_parameters"] == "cluster_analysis":
        lineage = "".join("".join(path.split("/")[-2]).split("-")[-1])
        busco_lineage(path, lineage, specific_name, db_path, out_dir, snakemake.threads, snakemake.log, copy=True)

    elif "auto-lineage-prok" in snakemake.config["busco_parameters"]:
        busco_autolineage(path, specific_name, db_path, out_dir, snakemake.threads, snakemake.log, copy=True)

    else:
        lineage = snakemake.config["busco_parameters"]
        busco_lineage(path, lineage, specific_name, db_path, out_dir, snakemake.threads, snakemake.log, copy=True)

else:
    logger.warning("Choose valid option for save_busco_output!")
    exit(1)