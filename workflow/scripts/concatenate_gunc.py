import pandas as pd
import glob

GUNC_output = pd.DataFrame()
for i in snakemake.input:
    tmp_df = pd.read_csv(i, sep="\t")
    GUNC_output = pd.concat([GUNC_output, tmp_df], ignore_index=True)
    GUNC_output.to_csv("results/GUNC_output.tsv", sep="\t")
