import os
import pandas as pd


BUSCO_output = pd.DataFrame()
for i in snakemake.input:
    command = f"grep -v '.snakemake' {i} > tmpfile && mv tmpfile {i}"
    os.system(command)
    tmp_df = (
        pd.read_csv(i, sep="\t")
        .reset_index()
        .shift(periods=3, axis=1)
        .dropna(axis=1)
    )
    BUSCO_output = pd.concat([BUSCO_output, tmp_df], ignore_index=True)
BUSCO_output.eval(
    "Completeness = (Complete + 0.5 * Fragmented) / n_markers", inplace=True
)
BUSCO_output.eval("Contamination = Duplicated / n_markers", inplace=True)
BUSCO_output.eval(
    "Quality_score = Completeness - 5*Contamination - (Fragmented/n_markers)",
    inplace=True,
)
BUSCO_output.to_csv("results/BUSCO_output.tsv", sep="\t")
