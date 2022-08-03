import glob
import pandas as pd
from snakemake.shell import shell
import os
import logging, traceback
try:
    from tqdm import tqdm
except ImportError:
    tqdm = []

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

# function to split list into chunks of size n
def chunks(lst, n):
    for i in range(0, len(lst), n):
        yield lst[i : i + n]


# get list of all genomes in directory specified in default_config
if snakemake.config['input_text_file'] == "":
    argument1 = f"{snakemake.config['genomes']}/*{snakemake.config['format1']}"
    argument2 = f"{snakemake.config['genomes']}/*{snakemake.config['format2']}"
    glob_argument1 = set(glob.glob(argument1))
    glob_argument2 = set(glob.glob(argument2))
    all_genomes_in_dir = list(glob_argument1.union(glob_argument2))
else:
    all_genomes_in_dir = []
    with open(snakemake.config['input_text_file'],'r') as f:
        lines = f.readlines()
        for line in lines:
            path = line.strip("\n")
            all_genomes_in_dir.append(path)
batch_size = int(snakemake.config["batch_size"])
N_batches = len(all_genomes_in_dir) // batch_size + 1
os.makedirs("quality_filtering/intermediate_results/genome_list/", exist_ok=True)
# if no other table is specified in default_config
if snakemake.config["busco_parameters"] != "cluster_analysis":
    logger.info(
        "No previous BUSCO output found! Proceeding with creating random batches . . ."
    )
    counter = 0
    generator = chunks(all_genomes_in_dir, batch_size)
    for batch in generator:
        key = snakemake.config['busco_parameters']
        counter += 1
        with open(
            f"quality_filtering/intermediate_results/genome_list/batch_number_{counter}{key}.txt",
            "w",
        ) as f:
            f.write("\n".join(batch))
else:
    logger.info(
        "Found previous BUSCO_output.tsv and table with information about clusters! Using it to create lineage-specific batches . . ."
    )
    all_genomes = pd.read_csv(
        snakemake.config["all_genomes"], sep="\t", index_col=0
    ).reset_index()
    BUSCO_output = pd.read_csv(snakemake.config["busco_output"], sep="\t", index_col=0)
    list_of_names_all_genomes = []
    os.makedirs("quality_filtering/intermediate_results/by_cluster/", exist_ok=True)
    # rename files to be more universal
    for i in range(len(all_genomes)):
        new_name = "_".join(
            "".join(
                all_genomes.iloc[i]["genome_file"].rstrip(".gz").split(".")[:-1]
            ).split("_")[:2]
        )
        list_of_names_all_genomes.append(new_name)
    all_genomes["genome_file_base"] = list_of_names_all_genomes
    list_of_names_BUSCO_output = []
    for n in range(len(BUSCO_output)):
        new_name = "_".join(
            BUSCO_output.iloc[n]["sample_name"].strip("[").strip("]").split("_")[:2]
        )
        list_of_names_BUSCO_output.append(new_name)
    BUSCO_output["genome_file_base"] = list_of_names_BUSCO_output
    # merge to create table that links cluster information with every genome in all_genomes
    merged = BUSCO_output.merge(all_genomes, on="genome_file_base")
    cluster_info = snakemake.config['cluster_info']
    big_table = all_genomes.merge(merged, on=cluster_info, how="left")
    # iterate over every cluster and create .txt files where one lineage = one .txt file
    # for i in range(1, len(big_table[cluster_info].unique())):
    #     df = big_table.query("@cluster_info == @i")
    for i, (_ , df) in tqdm(enumerate(big_table.groupby(cluster_info))):
        list_of_samples = []
        for sample_name in df["genome_file_base_x"]:
            list_of_samples.append(sample_name)
        try:
            if type(df.iloc[0]["lineage"]) == str:
                txt_folder = "quality_filtering/intermediate_results/by_cluster/"
                os.makedirs(txt_folder, exist_ok=True)
                txt_file = df.iloc[0]["lineage"]
                with open(f"{txt_folder}/links_{txt_file}.txt", "a") as f:
                    f.write("\n".join(list_of_samples) + "\n")
        except:
            pass
    logger.info(
        "Finished creating .txt files based on table with all genomes! Proceeding with loading input files . . ."
    )
    # Read those .txt files, create dictionary and that dataframe where columns = lineages and rows = samples
    dictionary_of_lineages = {}
    # glob_argument = "../../../HumGut/HumGut_all/by_cluster/*.txt"
    glob_argument = "quality_filtering/intermediate_results/by_cluster/*.txt"
    dictionary = {}
    for path in glob.glob(glob_argument):
        list_of_genomes = []
        lineage = "_".join("".join(path.rstrip(".txt").split("/")[-1]).split("_")[1:])
        with open(path, "r") as path_file:
            genomes = path_file.readlines()
            for genome in genomes:
                genome_stripped = genome.strip("\n")
                list_of_genomes.append(genome_stripped)
        dictionary[lineage] = list_of_genomes
    df_from_dict = pd.DataFrame.from_dict(dictionary, orient="index").T
    list_of_lineages = list(df_from_dict)
    for lineage in list_of_lineages:
        dictionary_of_lineages[lineage] = []
    # go through every sample in config['genomes'] dir and match the name with name in df_from_dict.
    # If there is a match, get name of the column (lineage)
    for sample in all_genomes_in_dir:
        sample_name = "_".join(
            "".join("".join(sample.rstrip(".gz").split("/")[-1]).split(".")[:-1]).split(
                "_"
            )[:2]
        ).strip("\n")
        try:
            lineage_of_sample = df_from_dict.columns[
                df_from_dict.isin([sample_name]).any()
            ][0]
            dictionary_of_lineages[lineage_of_sample].append(sample)
        except:
            logger.info(f"Cannot find lineage information for sample {sample_name}! Proceeding without it ...")
            pass
    counter = 0
    for key in list(dictionary_of_lineages.keys()):
        generator = chunks(dictionary_of_lineages[key], batch_size)
        for batch in generator:
            counter += 1
            # os.system(f"touch genome_list/batch_number_{counter}_lineage_{key}.txt")
            with open(
                f"quality_filtering/intermediate_results/genome_list/batch_number_{counter}-{key}.txt",
                "w",
            ) as f:
                f.write("\n".join(batch))
    logger.info(f"Finished creating batches!")
