from glob import glob
import os, sys
import logging, traceback
from itertools import repeat
from snakemake.shell import shell
import subprocess
import collections

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
    from tqdm import tqdm
except ImportError:
    tqdm= list

    logging.warning("optional package tqdm is not installed")



def shell(cmd):
    """ My shell command that reports error to the logger
    """
    sp= subprocess.run(cmd,shell=True,capture_output=True,text=True)

    #if sp.returncode!=0:
    #    raise sp.CalledProcessError(sp.retcode, cmd,stderr=sp.stderr,stdout=sp.sdtout)


def splitext_ignore_gz(path):
    basename,ext = os.path.splitext(path)

    if ext=='.gz':
        basename,ext = os.path.splitext(basename)

    return basename, ext




def any2fasta(path_list, intermediate_fasta_dir):
    global sample_names
#read paths from input file
    with open(path_list, "r") as path_file:
        paths= path_file.read().split()

    sample_names=[]

    for file in tqdm(paths):
        filename, extension= splitext_ignore_gz(file)
        sample_name = os.path.basename(filename)

        sample_names.append(sample_name)

        if (extension==".gff") or file.endswith(".gz"):
            logging.info(
                f"Sample {sample_name} is eather zipped or not in fasta format! Proceeding with any2fasta. . ."
            )
            shell(
                f"any2fasta {file} > {intermediate_fasta_dir}/{sample_name}.fasta"
            )
        else:
            logging.info(
                f"Sample {sample_name} is in fasta format! Proceeding with creating symlink. . ."
            )

            os.symlink(file, f"{intermediate_fasta_dir}/{sample_name}.fasta" )

### Run prodigal


from multiprocessing import Pool


def run_prodigal(input_fasta, parameters, out_dir, prodigal_output, dictionary):#snakemake.params.parameters ):
    """This rule has only one input argument othersise the multithreding becomes a bit more complicated
    """
    # run prodigal on one sample
    # stderr get raised and saved to log file
    sample_name = os.path.splitext(os.path.basename(input_fasta))[0]
    # try:
    path = dictionary[sample_name]
    os.makedirs(f"{prodigal_output}/faa/{path}", exist_ok=True)
    os.makedirs(f"{prodigal_output}/fna/{path}", exist_ok=True)
    logging.info(f"Run prodigal on sample {sample_name}")
    shell(f"prodigal {parameters} -i {input_fasta} -a {prodigal_output}/faa/{path}/{sample_name}.faa  -d {prodigal_output}/fna/{path}/{sample_name}.fna > /dev/null")
    # except:
    #     logging.info(f"Couldn't run prodigal on {sample_name}. Countuing without it ...")

def run_multiple_prodigal(input_dir, threads, extension=".fasta"):

    # read input list
    genomes_fastas = glob(os.path.join(input_dir, "*" + extension))
    assert len(genomes_fastas) > 0


    pool = Pool(int(threads))
    pool.starmap(run_prodigal, zip(genomes_fastas, repeat(snakemake.params), repeat(out_dir), repeat(prodigal_output), repeat(dictionary)))

def prepare_faa_paths(paths):
    dictionary = collections.defaultdict()
    for path in paths:
        full_sample = os.path.basename(path)
        sample = full_sample.rstrip(snakemake.config['format1']).rstrip(snakemake.config['format2'])
        path_for_faa = "/".join(path.split("/")[-4:-1]).rstrip(snakemake.config['format1']).rstrip(snakemake.config['format2'])
        dictionary[sample] = path_for_faa
    
    return dictionary

if __name__ == "__main__":

    logging.basicConfig(
        filename=snakemake.log[0],
        level=logging.INFO,
        format="%(asctime)s %(message)s",
        datefmt="%Y-%m-%d %H:%M:%S",
    )
    out_dir = "quality_filtering/"
    # create fasta dir
    intermediate_fasta_dir = f"{snakemake.config['temporary_dir']}/intermediate_results/fasta/{snakemake.wildcards.counter}-{snakemake.wildcards.lineage}"
    os.makedirs(intermediate_fasta_dir, exist_ok=True)

    prodigal_output = snakemake.config.get('directory_faa', 'quality_filtering/faa_files')
    for subfolders in ["faa","fna"]:
        os.makedirs(f"{prodigal_output}/{subfolders}",exist_ok=True)

    if snakemake.config['input_text_file'] != "":
        with open(snakemake.config['input_text_file'],'r') as f:
            lines = f.readlines()
            dictionary = prepare_faa_paths(lines)
    elif snakemake.config['genomes'] != "":
        paths = glob(f"{snakemake.config['genomes']}/*{snakemake.config['format1']}")
        paths.extend(glob(f"{snakemake.config['genomes']}/*{snakemake.config['format2']}"))
        dictionary = prepare_faa_paths(paths)

    any2fasta(snakemake.input[0],intermediate_fasta_dir)
    run_multiple_prodigal(intermediate_fasta_dir, threads = snakemake.threads)



    ### Create symlinks for faa file

    # get input_dir relative from output dir to create a relative symlink
    os.makedirs(snakemake.output[0],exist_ok=True)
    input_dir_rel = os.path.relpath(f"{prodigal_output}/faa/", snakemake.output[0])
    for sample_name in sample_names:
        part_path = dictionary[sample_name]
        os.symlink(f"{input_dir_rel}/{part_path}/{sample_name}.faa",
                        f"{snakemake.output[0]}/{sample_name}.faa" )
