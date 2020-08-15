import os
from snakemake.io import expand

def read_samples(samples_file, data_dir):
    sample_list = [ line.strip() for line in open(samples_file, 'rt') ]
    sample_list = [ line for line in sample_list if line ]   # remove empty lines
    # verify that all genome files exist -
    for filename in sample_list:
        fullpath = os.path.join(data_dir, filename)
        if not os.path.exists(fullpath):
            print(f'** ERROR: genome file {filename} does not exist in {data_dir}')
            if strict_mode:
                print('** exiting.')
                sys.exit(-1)
    return sample_list

def check_params(config):
    pass
    # to do: check that scaled, ksizes, etc make sense for sourmash / for each alpha
    #ksize = config['ksize']
    #try:
    #    ksize = int(ksize)
    #    if ksize < 15 or ksize > 101:
    #        raise ValueError
    #except ValueError:
    #    print('** ERROR: ksize should be a number between 15 and 101.')
    #    print('** (it must also match the query database ksize value)')
    #    if strict_mode:
    #        sys.exit(-1)


def generate_targets(config, samples, output_dir=""):
    # get pipeline we're using (default = taxonomic_classification_gtdb)
    pipeline = config["pipeline"]
    # find required db's for this pipeline
    databases = config["pipelines"][pipeline].get("databases", [])

    # Databases:: generate targets for each database
    database_dir = config["database_dir"] # to do: make/use a sanitize path function to remove `~` and get abspath here
    database_targets=[]
    for db in databases:
        ## to do  -- cleaner/clearer/better db specification in yaml file! Not so much nesting?
        ## instead of nesting protein - name them uniquely? handier for matching exact db in config.yaml
        database_targets+= [os.path.join(database_dir, config["databases"][db]["protein"]["k11"]["sbt"])]

    # Pipeline:: find steps in this pipeline
    input_type = config["input_type"]
    steps = []
    # if nucleotide input, run both protein and nucl steps, else just run protein steps
    if input_type in ["protein", "nucleotide"]:
        if input_type == "nucleotide":
            steps  = config["pipelines"][pipeline]["steps"]["nucleotide"]
        steps += config["pipelines"][pipeline]["steps"]["protein"]
    else:
        print(f'** ERROR: input type {input_type} must be either "protein" or "nucleotide"')
        if strict_mode:
            print('** exiting.')
            sys.exit(-1)

    # generate targets for each step
    pipeline_targets=[]
    for step in steps:
        step_outdir = config[step]["output_dir"]
        step_files = config[step]["output_files"]
        step_databases = config[step].get("databases", [])

        # fill variables in the output filenames
        for stepF in step_files:
            pipeline_targets += expand(os.path.join(output_dir, step_outdir, stepF), sample=samples, database=step_databases)

    # optionally only return databases or pipeline to enable running these as separate bits of the workflow
    targets = database_targets + pipeline_targets
    return targets


# for nucleotide vs protein: could specify steps as nucleotide or protein, choose steps based on config["input_type"]
