import os
import sys
import pandas as pd
from snakemake.io import expand


def find_input_file(filename, name="input", add_paths=[], add_suffixes = ['.yaml', '.yml'], verbose = False):
    # for any file specified via command line, check if it exists at the current path, if not, try some other paths before returning  a helpful error
    found_file = None
    filename = sanitize_path(filename) # handle ~!
    paths_to_try = ['', os.getcwd(), os.path.dirname(os.path.abspath(__file__)), os.path.dirname(os.path.dirname(os.path.abspath(__file__)))] + add_paths
    suffixes_to_try = [''] + add_suffixes

    if os.path.exists(filename) and not os.path.isdir(filename):
        found_file = os.path.realpath(filename)
    else:
        for p in paths_to_try:
            for s in suffixes_to_try:
                tryfile = os.path.join(p, filename+ s)
                if os.path.exists(tryfile) and not os.path.isdir(tryfile):
                    found_file = os.path.realpath(tryfile)
                    break
    assert found_file, f'Error, cannot find specified {name} file {filename}\n\n\n'
    if verbose:
        sys.stderr.write(f'\tFound {name} at {found_file}\n')
    return found_file

def read_samples(samples_file, data_dir, strict_mode=False):
    samples_file = find_input_file(samples_file)
    if '.tsv' in samples_file or '.csv' in samples_file:
        separator = '\t'
        if '.csv' in samples_file:
            separator = ','
        try:
            samples = pd.read_csv(samples_file, dtype=str, sep=separator, names = ["sample", "filename"])
        except Exception as e:
            sys.stderr.write(f"\n\tError: {samples_file} file is not properly formatted. Please fix.\n\n")
            print(e)
    elif '.xls' in samples_file:
        try:
            samples = pd.read_excel(samples_file, dtype=str, names = ["sample", "filename"])
        except Exception as e:
            sys.stderr.write(f"\n\tError: {samples_file} file is not properly formatted. Please fix.\n\n")
            print(e)
    else:
        sample_list = [ line.strip() for line in open(samples_file, 'rt') ]
        sample_list = [ line for line in sample_list if line ]   # remove empty lines
        samplesD = {"sample":sample_list,"filename":sample_list} # maybe later try removing *fa.gz or the like
        samples = pd.DataFrame(samplesD)

    samples.set_index("sample", inplace=True)

    # Now, verify that all genome files exist
    data_dir = sanitize_path(data_dir)
    sample_list = samples["filename"].tolist()
    for filename in sample_list:
        fullpath = os.path.join(data_dir, filename)
        if not os.path.exists(fullpath):
            print(f'** ERROR: genome file {filename} does not exist in {data_dir}')
            if strict_mode:
                print('** exiting.')
                sys.exit(-1)

    return samples


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


def check_input_type(config, strict_mode=False):
    """
    Check the input type provided in the configfile
    """
    alphabet_info = config["alphabet_info"]
    input_type = config.get("input_type")
    if input_type:
        if input_type in ["protein", "nucleotide"]:
            if input_type == "protein" and "nucleotide" in alphabet_info.keys():
                print(f'** input type is protein. Turning off nucleotide searching.')
                del alphabet_info["nucleotide"]
                config["alphabet_info"] = alphabet_info
        else:
            print(f'** ERROR: input type {input_type} must be either "protein" or "nucleotide"')
            if strict_mode:
                print('** exiting.')
                sys.exit(-1)
    else:
        print(f'** ERROR: "input_type: protein" or "input_type: nucleotide" must be specified in the config file')
        print('** exiting.')
        sys.exit(-1)
    return config

def check_and_set_alphabets(config, strict_mode=False):
    """
    Choose the alphabets for sketching and search
    """
    alphaInfo={}

    default_alphabets = config["alphabet_defaults"]
    alphabets = config.get("alphabets", default_alphabets)
    # check each alphabet
    for alpha in alphabets:
        if alpha not in default_alphabets.keys():
            print(f'** ERROR: alphabet {alpha} is invalid. Valid alphabets are: {", ".join(default_alphabets)}')
            if strict_mode:
                print('** exiting.')
                sys.exit(-1)
            else:
                print(f'Strict mode is off: attempting to continue. Removing {alpha} from alphabet list.')
                alphabets.remove(alpha)
        else:
            # first, set defaults
            alphaInfo[alpha] = default_alphabets[alpha]
            # if override values are provided in the config, replace defaults
            if isinstance(alphabets, dict):
                alphaInfo[alpha]["ksizes"] = alphabets[alpha].get("ksizes", default_alphabets[alpha]["ksizes"])
                alphaInfo[alpha]["scaled"] = alphabets[alpha].get("scaled", default_alphabets[alpha]["scaled"])
    if alphaInfo:
        config["alphabet_info"] = alphaInfo
        config = check_input_type(config)
    else:
        print(f'** ERROR: no valid alphabets remain')
        print('** exiting.')
        sys.exit(-1)

    return config



def sanitize_path(path):
    # expand `~`, get absolute path
    path = os.path.expanduser(path)
    path = os.path.abspath(path)
    return path



def check_dbinfo_exists(db_name, all_dbinfo, strict_mode=False):
    dbinfo_exists = False
    available_databases = all_dbinfo.keys()
    if db_name in available_databases:
        # maybe check alphas and ksizes here?
        dbinfo_exists=True
    else:
        print(f'** ERROR: database {db_name} does not have any associated info. Available databases are: {", ".join(available_databases)}')
        if strict_mode:
            print('** exiting.')
            sys.exit(-1)
        else:
            print('Strict mode is off: attempting to continue. Removing {db_name} from search databases list.')
    return dbinfo_exists


def check_user_databases_and_set_info(config, strict_mode=False):
    databases=[]

    # get information for available databases (ksizes, file dl info, etc)
    db_info = config["database_info"]

    # find the default databases for this pipeline
    pipeline = config["pipeline"]
    default_databases = config["pipelines"][pipeline].get("databases", [])

    # use pipeline dbs unless user turns them off
    no_use_defaults = config.get("turn_off_default_databases", False)
    if not no_use_defaults:
        # check that we have the info for these default databases
        for db in default_databases:
            if check_dbinfo_exists(db, db_info):
                databases.append(db)

    ## add user databases if info is available
    user_dbs = config.get("search_databases", [])

    for db in user_dbs:
        if check_dbinfo_exists(db, db_info):
            # add database to list of databases to search
            databases.append(db)

    if not databases:
        print(f'** ERROR: no valid databases selected. Please choose from the available databases or do not disable default databases.')
        print('** exiting.')
        sys.exit(-1)

    ## TO DO: ##
    # check for alpha-ksize for each db; only run alpha-ksizes that exist in dbs

    config["databases"] = databases

    # sanitize database directory path
    config["database_dir"] = sanitize_path(config["database_dir"])

    return config


def integrate_user_config(config):
    config = check_and_set_alphabets(config)
    pipeline = config["pipeline"]
    if config["pipelines"][pipeline]["databases_required"]:
        config["get_databases"] = True
        config = check_user_databases_and_set_info(config)
    else:
        config["get_databases"] = False
    return config



def generate_database_targets(config, also_return_database_names=False):
    database_targets, database_names=[],[]
    ## integrate user settings and inputs
    config = integrate_user_config(config)
    ## What alphabets are we using? ##
    alphabet_info = config["alphabet_info"]
    ## What databases are we using? ##
    databases = config["databases"]
    # get all database details
    database_info = config["database_info"]

    # default filenaming for each database
    # variables: db_name, alphabet, ksize, db_type, suffix
    db_target_templates= config["database_target_template"]
    info_templates = db_target_templates["info_csv"]
    db_templates = db_target_templates["database"]

    # iterate through dbinfo and build targets for the alphabets we're using
    for db in databases:
        db_targs,db_names=[],[]
        db_info = config["database_info"][db]
        for db_alphabet, db_alpha_info in db_info["alphabets"].items():
            if db_alphabet in alphabet_info.keys():
                for db_ksize, dbs in db_alpha_info.items():
                    ksize_int = int(db_ksize[1:]) # db_ksize is string w/format: k{ksize}
                    # only build target if db has matching ksize available.
                    # todo: also handle scaled here???
                    if ksize_int in alphabet_info[db_alphabet]["ksizes"]:

                        for db_type in dbs.keys():
                            suffix = config["database_suffixes"][db_type]
                            db_filenames = expand(db_templates, db_name=db,alphabet=db_alphabet, ksize=db_ksize, db_type=db_type, suffix=suffix)
                            # generate db_name, needed for workflow targets. sigh, don't like this - do it better.
                            end = f".{db_type}.{suffix}"
                            names = [fn.rsplit(end)[0] for fn in db_filenames]
                            db_targs+=db_filenames
                            db_names+=names
        # if we have any targets for this database name, also grab the info csv target
        if db_targs:
            # also get the db info csv
            db_info = expand(info_templates, db_name=db)
            db_targs+=db_info
        # add targets for this database
        database_targets+=db_targs
        database_names+=db_names

    database_dir=config["database_dir"]
    final_db_targs = [os.path.join(database_dir, x) for x in database_targets]
    if also_return_database_names:
        return final_db_targs, database_names

    return final_db_targs


def generate_targets(config, samples, output_dir="", generate_db_targets=False):
    pipeline_targets=[]
    ## integrate user settings and inputs
    config = integrate_user_config(config)
    ## What alphabets are we using? ##
    alphabet_info = config["alphabet_info"]
    ## set run basename
    basename = config.get("basename", "thumper-output")

    ## What databases are we using? ##
    # Pipeline:: find steps in this pipeline
    pipeline= config["pipeline"]
    db_required = config["pipelines"][pipeline]["databases_required"]
    if db_required:
        database_targets, database_names = generate_database_targets(config, also_return_database_names=True)
    else:
        generate_db_targets=False
        database_targets,database_names=[],[]
    index_names = []
    if pipeline == "generate_index":
        for alpha, alphaInfo in alphabet_info.items():
            index_names+=expand("{basename}.{alpha}-k{ksize}.scaled{scaled}", basename=basename, alpha=alpha, ksize=alphaInfo["ksizes"], scaled=alphaInfo["scaled"])

    # generate targets for each step
    steps = config["pipelines"][pipeline]["steps"]
    for step in steps:
        step_outdir = config[step]["output_dir"]
        step_files = config[step]["output_files"]

        # fill variables in the output filenames
        for stepF in step_files:
            pipeline_targets += expand(os.path.join(output_dir, step_outdir, stepF), sample=samples, database=database_names, basename=basename, db_name=config.get("databases", []), index=index_names)

    if generate_db_targets:
        targets = database_targets + pipeline_targets
        return targets

    return pipeline_targets

