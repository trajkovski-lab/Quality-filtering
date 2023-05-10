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

# Script starts here

tables = [pd.read_csv(f"{input_dir}/batch_summary.txt", sep='\t', on_bad_lines='warn') for input_dir in snakemake.input]
df = pd.concat(tables, ignore_index=True)

df = df.drop(['Scores_archaea_odb10','Scores_bacteria_odb10','Scores_eukaryota_odb10'], axis=1)

df = df.eval("Completeness = (Complete + 0.5 * Fragmented)/100")
df = df.eval("Contamination = Duplicated / 100")
df = df.eval("Quality_score = Completeness - 5*Contamination - (Fragmented/100)")

df.to_csv(snakemake.output[0], sep="\t", index=False)