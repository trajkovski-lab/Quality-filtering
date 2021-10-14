import glob

def chunks(lst, n):
    for i in range(0, len(lst), n):
        yield lst[i:i + n]


argument1 = f"{snakemake.config['genomes']}/*{snakemake.config['format1']}"
argument2 = f"{snakemake.config['genomes']}/*{snakemake.config['format2']}"
glob_argument1 = set(glob.glob(argument1))
glob_argument2 = set(glob.glob(argument2))
all_genomes = list(glob_argument1.union(glob_argument2))
batch_size = snakemake.config['batch_size']
generator = chunks(all_genomes, batch_size)
N_batches = len(all_genomes) // batch_size + 1
i=0
while i < N_batches:
    for batch in generator:
        with open(snakemake.output[i], 'w') as f:
            f.write("\n".join(batch))
            i += 1
'''
n = len(snakemake.output)
step = len(all_genomes) // n
for i in range(n):
    with open(snakemake.output[i], "w") as f:
        f.write("\n".join(all_genomes[(i * step) : ((i + 1) * step)]))
        if i <= (len(all_genomes) % n):
            f.write("\n" + all_genomes[-i])
'''
