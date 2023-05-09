import os
import sys
import pandas as pd
import glob
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

final_df = pd.DataFrame()
for input_dir in snakemake.input:
    table = f"{input_dir}/quality_report.tsv"
    logger.info(f"Found {table} table! Concatenating it to checkm2_output.tsv . . .")
    checkm2_intermediary = pd.read_csv(table, sep="\t")
    final_df = final_df.append(checkm2_intermediary, ignore_index=True)
    final_df.to_csv(snakemake.output[0], sep="\t", index=False)
