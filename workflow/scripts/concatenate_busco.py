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
for table in snakemake.input:
    logger.info(f"Found {table} table! Concatenating it to BUSCO_output.tsv . . .")
    busco_intermediary = pd.read_csv(table, sep="\t", index_col=0)
    final_df = final_df.append(busco_intermediary, ignore_index=True)
    final_df.to_csv(snakemake.output[0], sep="\t")
