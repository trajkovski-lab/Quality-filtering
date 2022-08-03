import logging, traceback
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

## Start script
import os

# remove tailing '/' if it exitst. omitts snakemake to give warnings
config["database_folder"] = config["database_folder"].strip("/")

def calculate_N_batches():

    if 'N_batches' not in config:

        logger.info("Calculate number of batches")

        if config['input_text_file'] == "":
            files = glob.glob(os.path.join( config['genomes'] ,"*" + config["format1"]))

            if config["format2"] is not None:
                files.extend( glob.glob(os.path.join( config['genomes'] ,"*" + config["format2"])) )
        
            N_files = len(set(files))
        else: 
            N_files = 0
            with open(config['input_text_file'],'r') as f:
                lines = f.readlines()
                for line in lines:
                    genome_path = line.rstrip('\n')
                    if os.path.exists(genome_path) and os.path.getsize(genome_path) > 0:
                        N_files +=1
                    else:
                        logger.warning(f"Genome {genome_path} doesn't exist or is empty!")

        if N_files==0 :
                logger.warning("Didn't found genomes in folder {genomes} with extension {format1} or {fomrat2}".format(**config))
                exit(1)
        
        config['N_batches'] = N_files // config['batch_size'] + 1
        logger.info("Process {N_files} genomes in {N_batches} batches of size <= {batch_size}".format(N_files= N_files, **config))
    return config['N_batches']


scattergather:
    split= calculate_N_batches()


configfile: f"{os.path.dirname(workflow.snakefile)}/../config/default_config.yaml"
