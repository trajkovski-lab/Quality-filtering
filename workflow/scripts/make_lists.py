import glob

glob_argument = f"{snakemake.config['genomes']}/*{snakemake.config['format']}"
all_genomes = glob.glob(glob_argument)
n = len(snakemake.output)
step = len(all_genomes) // n
for i in range(n):
    with open(snakemake.output[i], "w") as f:
        f.write(" ".join(all_genomes[(i * step) : ((i + 1) * step)]))
        if i <= (len(all_genomes) % n):
            f.write(" " + all_genomes[-i])
