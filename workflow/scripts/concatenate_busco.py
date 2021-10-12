import os
import sys
import pandas as pd
import glob
final_df = pd.DataFrame()
temporary_df = pd.DataFrame()
with open(snakemake.log[0], "w") as f:
    sys.stderr = sys.stdout = f
    print("odmah ispod")
    for i in snakemake.input:
        glob_argument = f"{i}/*/short_summary.specific.*.*.txt"
        print(i)
        for path in glob.glob(glob_argument):
            sample_name = path.split("/")[-1].split(".")[-2]
            lineage = path.split("/")[-1].split(".")[-3]
            with open(path, 'r') as f:
                lines = f.readlines()
                counter = -1
                for line in lines:
                    counter += 1
                    if counter == 8:
                        numbers = [int(word) for word in line.lstrip("\t").split() if word.isdigit()]
                        complete_busco = int(*numbers)
                    if counter == 9:
                        numbers = [int(word) for word in line.lstrip("\t").split() if word.isdigit()]
                        single_busco = int(*numbers)
                    if counter == 10:
                        numbers = [int(word) for word in line.lstrip("\t").split() if word.isdigit()]
                        duplicated_busco = int(*numbers)
                    if counter == 11:
                        numbers = [int(word) for word in line.lstrip("\t").split() if word.isdigit()]
                        fragmented_busco = int(*numbers)
                    if counter == 12:
                        numbers = [int(word) for word in line.lstrip("\t").split() if word.isdigit()]
                        missing_busco = int(*numbers)
                    if counter == 13:
                        numbers = [int(word) for word in line.lstrip("\t").split() if word.isdigit()]
                        total_busco = int(*numbers)
                        series = pd.Series([sample_name,lineage,complete_busco, single_busco, duplicated_busco,fragmented_busco,missing_busco,total_busco], index = ["sample_name","lineage","complete_busco","single_busco","duplicated_busco","fragmented_busco","missing_busco","total_busco"])
            temporary_df = temporary_df.append(series, ignore_index= True)
        final_df = final_df.append(temporary_df)
    temporary_df.eval("Completeness = (complete_busco + 0.5 * fragmented_busco) / total_busco", inplace=True)
    temporary_df.eval("Contamination = duplicated_busco / total_busco", inplace=True)
    temporary_df.eval("Quality_score = Completeness - 5*Contamination - (fragmented_busco/total_busco)",inplace=True)
    temporary_df.to_csv("results/BUSCO_output.tsv", sep = '\t')
