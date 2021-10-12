

# remove tailing '/' if it exitst. omitts snakemake to give warnings
config["database_folder"] = config["database_folder"].strip("/")

def calculate_N_batches():

    if 'N_batches' not in config:

        logger.info("Calculate number of batches")

        N_files = len( glob.glob(os.path.join( config['genomes'] ,"*" + config["format1"])) )

        if config["format2"] is not None:
            N_files += len( glob.glob(os.path.join( config['genomes'] ,"*" + config["format2"])) )


        if N_files==0 :
            logger.critical("Didn't found genomes in folder {genomes} with extension {format1} or {fomrat2}".format(**config))
            exit(1)

        config['N_batches'] = N_files // config['batch_size'] + 1

        logger.info("Process {N_files} genomes in {N_batches} of size <= {batch_size}".format(N_files= N_files, **config))

    return config['N_batches']


scattergather:
    split= calculate_N_batches()


configfile: f"{os.path.dirname(workflow.snakefile)}/../config/default_config.yaml"
