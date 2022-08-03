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
    table = f"{input_dir}/batch_summary.txt"
    logger.info(f"Found {table} table! Concatenating it to BUSCO_output.tsv . . .")
    busco_intermediary = pd.read_csv(table, sep="\t").dropna().drop(['Scores_archaea_odb10','Scores_bacteria_odb10','Scores_eukaryota_odb10'],axis=1)
    busco_intermediary.eval("Completeness = (Complete + 0.5 * Fragmented)/100",inplace=True)
    busco_intermediary.eval("Contamination = Duplicated / 100",inplace=True)
    busco_intermediary.eval("Quality_score = Completeness - 5*Contamination - (Fragmented/100)",inplace=True)
    final_df = final_df.append(busco_intermediary, ignore_index=True)
    final_df.to_csv(snakemake.output[0], sep="\t", index=False)
