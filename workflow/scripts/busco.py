import glob
import os
import pandas as pd
for i in snakemake.input.fasta_dir:
    temporary_df = pd.DataFrame()
    glob_argument = f"{i}/*.faa"
    for path in glob.glob(glob_argument):
        sample_name = "".join(path.split("/")[-1].split(".")[:-1])
        command = f"busco -i {path} {snakemake.config['busco_parameters']} -m protein -o {sample_name} --out_path {snakemake.params.working_dir} -c {snakemake.threads} &> {snakemake.log}"
        os.system(command)
        glob_result = f"{snakemake.params.working_dir}/*/short_summary.specific.*.*.txt"
    for path_result in glob.glob(glob_result):
        sample_name = path_result.split("/")[-1].split(".")[-2]
        lineage = path_result.split("/")[-1].split(".")[-3]
        with open(path_result, 'r') as f:
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
        temporary_df.eval("Completeness = (complete_busco + 0.5 * fragmented_busco) / total_busco", inplace=True)
        temporary_df.eval("Contamination = duplicated_busco / total_busco", inplace=True)
        temporary_df.eval("Quality_score = Completeness - 5*Contamination - (fragmented_busco/total_busco)",inplace=True)
    file_name = f"results/busco/{snakemake.wildcards.scatteritem}.tsv"
    temporary_df.to_csv(file_name, sep = '\t')
