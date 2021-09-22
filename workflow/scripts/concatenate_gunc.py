import pandas as pd
import glob

GUNC_output = pd.DataFrame()
for i in snakemake.input:
    glob_argument = f"{i}/*.tsv"
    for path_to_tsv in glob.glob(glob_argument):
        tmp_df = pd.read_csv(path_to_tsv, sep="\t")
        GUNC_output = pd.concat([GUNC_output, tmp_df], ignore_index=True)
        GUNC_output.to_csv("results/GUNC_output.tsv", sep="\t")
