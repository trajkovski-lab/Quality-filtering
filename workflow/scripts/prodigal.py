from glob import glob
#from snakemake.shell import shell
import os, sys
import logging, traceback

import subprocess

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


try:
    import tqdm
except ImportError:
    tqdm= list

    logger.warning("optional package tqdm is not installed")



def shell(cmd):
    """ My shell command that reports error to the logger
    """
    sp= subprocess.run(command,shell=True,capture_output=True,text=True)

    if sp.returncode!=0:
        raise sp.CalledProcessError(sp.retcode, cmd,stderr=sp.stderr,stdout=sp.sdtout)


def splitext_ignore_gz(path):
    basename,ext = os.path.splitext(path)

    if ext=='.gz':
        basename,ext = os.path.splitext(basename)

    return basename, ext



# create fasta dir
intermediate_fasta_dir = f"{snakemake.config['temporary_dir']}/intermediate_results/fasta/{snakemake.wildcards.counter}-{snakemake.wildcards.lineage}"
os.makedirs(intermediate_fasta_dir)

# create faa dir
directory_faa = snakemake.config.get('directory_faa', 'quality_filtering/faa_files')
os.makedirs(directory_faa,exist_ok=True)




#read paths from input file
with open(snakemake.input[0], "r") as path_file:
    paths= path_file.read().split()



for file in tqdm(paths):


    filename, extension= splitext_ignore_gz(file)
    sample_name = os.path.basename(filename)


    if (extension==".gff") or file.endswith(".gz"):
        logger.info(
            f"Sample {sample_name} is eather zipped or not in fasta format! Proceeding with any2fasta. . ."
        )

        shell(
            f"any2fasta {file} > {intermediate_fasta_dir}/{sample_name}.fasta"
        )
    else:
        logger.info(
            f"Sample {sample_name} is in fasta format! Proceeding with creating symlink. . ."
        )

        os.path.symlink(os.path.abspath(file), f"{intermediate_fasta_dir}/{sample_name}.fasta" )


### Run prodigal


from multiprocessing import Pool


def run_prodigal(input_fasta,out_dir=directory_faa, parameters=snakemake.params.parameters ):
    """This rule has only one input argument othersise the multithreding becomes a bit more complicated
    """
    # run prodigal on one sample
    # stderr get raised and saved to log file
    sample_name = os.path.splitext(os.path.basename(input_fasta))[0]
    logging.info(f"Run prodigal on sample {sample_name}")

    shell(f"prodigal {parameters} -i {input_fasta} -a {out_dir}/{sample_name}.faa > /dev/null")



def run_multiple_prodigal(input_dir, threads=snakemake.threads, extension=".fasta"):

    # read input list
    genomes_fastas = glob(os.path.join(input_dir, "*" + extension))
    assert len(genomes_fastas) > 0


    pool = Pool(int(threads))
    pool.map(run_prodigal, genomes_fastas)


run_multiple_prodigal(intermediate_fasta_dir)
