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

GUNC_output = pd.DataFrame()
for i in snakemake.input:
    logger.info(f"Found {i} table! Concatenating it to GUNC_output.tsv . . .")
    tmp_df = pd.read_csv(i, sep="\t")
    GUNC_output = pd.concat([GUNC_output, tmp_df], ignore_index=True)
    GUNC_output.to_csv("quality_filtering/GUNC_output.tsv", sep="\t")
