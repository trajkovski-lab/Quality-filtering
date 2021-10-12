import glob

argument1 = f"{snakemake.config['genomes']}/*{snakemake.config['format1']}"
argument2 = f"{snakemake.config['genomes']}/*{snakemake.config['format2']}"
glob_argument1 = set(glob.glob(argument1))
glob_argument2 = set(glob.glob(argument2))
all_genomes = list(glob_argument1.union(glob_argument2))
n = len(snakemake.output)
step = len(all_genomes) // n
for i in range(n):
    with open(snakemake.output[i], "w") as f:
        f.write(" ".join(all_genomes[(i * step) : ((i + 1) * step)]))
        if i <= (len(all_genomes) % n):
            f.write(" " + all_genomes[-i])
