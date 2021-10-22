import glob
import os
import pandas as pd
from snakemake.shell import shell
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


for i in snakemake.input.fasta_dir:
    temporary_df = pd.DataFrame()
    glob_argument = f"{i}/*.faa"
    db_path = os.path.join(snakemake.config["database_folder"]
    for path in glob.glob(glob_argument):
        sample_name = ".".join(path.split("/")[-1].split(".")[:-1])
        if snakemake.config["busco_parameters"] == "cluster_analysis":
            lineage = "".join("".join(path.split("/")[-2]).split("-")[-1])
            shell(
                f"busco -i {path} -l {lineage} -m protein -o {sample_name} --download_path {db_path} --out_path {snakemake.params.working_dir} -c {snakemake.threads} &> {snakemake.log}"
            )
        elif "auto-lineage-prok" in snakemake.config["busco_parameters"]:
            shell(
                f"busco -i {path} --auto-lineage-prok -m protein -o {sample_name} --download_path {db_path} --out_path {snakemake.params.working_dir} -c {snakemake.threads} &> {snakemake.log}"
            )
        else:
            lineage = snakemake.config["busco_parameters"]
            shell(
                f"busco -i {path} -l {lineage} -m protein -o {sample_name} --download_path {db_path} --out_path {snakemake.params.working_dir} -c {snakemake.threads} &> {snakemake.log}"
            )
    glob_result = f"{snakemake.params.working_dir}/*/short_summary.specific.*.*.txt"
    for path_result in glob.glob(glob_result):
        sample_name = path_result.split("/")[-1].split(".")[-2]
        lineage = path_result.split("/")[-1].split(".")[-3]
        with open(path_result, "r") as f:
            lines = f.readlines()
            counter = -1
            for line in lines:
                counter += 1
                if counter == 8:
                    numbers = [
                        int(word)
                        for word in line.lstrip("\t").split()
                        if word.isdigit()
                    ]
                    complete_busco = int(*numbers)
                if counter == 9:
                    numbers = [
                        int(word)
                        for word in line.lstrip("\t").split()
                        if word.isdigit()
                    ]
                    single_busco = int(*numbers)
                if counter == 10:
                    numbers = [
                        int(word)
                        for word in line.lstrip("\t").split()
                        if word.isdigit()
                    ]
                    duplicated_busco = int(*numbers)
                if counter == 11:
                    numbers = [
                        int(word)
                        for word in line.lstrip("\t").split()
                        if word.isdigit()
                    ]
                    fragmented_busco = int(*numbers)
                if counter == 12:
                    numbers = [
                        int(word)
                        for word in line.lstrip("\t").split()
                        if word.isdigit()
                    ]
                    missing_busco = int(*numbers)
                if counter == 13:
                    numbers = [
                        int(word)
                        for word in line.lstrip("\t").split()
                        if word.isdigit()
                    ]
                    total_busco = int(*numbers)
                    series = pd.Series(
                        [
                            sample_name,
                            lineage,
                            complete_busco,
                            single_busco,
                            duplicated_busco,
                            fragmented_busco,
                            missing_busco,
                            total_busco,
                        ],
                        index=[
                            "sample_name",
                            "lineage",
                            "complete_busco",
                            "single_busco",
                            "duplicated_busco",
                            "fragmented_busco",
                            "missing_busco",
                            "total_busco",
                        ],
                    )
                    temporary_df = temporary_df.append(series, ignore_index=True)
        temporary_df.eval(
            "Completeness = (complete_busco + 0.5 * fragmented_busco) / total_busco",
            inplace=True,
        )
        temporary_df.eval(
            "Contamination = duplicated_busco / total_busco", inplace=True
        )
        temporary_df.eval(
            "Quality_score = Completeness - 5*Contamination - (fragmented_busco/total_busco)",
            inplace=True,
        )
    file_name = f"quality_filtering/intermediate_results/busco_tables/{snakemake.wildcards.counter}-{snakemake.wildcards.lineage}.tsv"
    temporary_df.to_csv(file_name, sep="\t")
