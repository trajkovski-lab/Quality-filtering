import os
import sys
import pandas as pd
import glob
final_df = pd.DataFrame()
for table in snakemake.input:
    busco_intermediary = pd.read_csv(table, sep = '\t', index_col = 0)
    final_df = final_df.append(busco_intermediary, ignore_index=True)
    final_df.to_csv(snakemake.output[0], sep = '\t')
