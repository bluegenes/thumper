import os
import sys
import pandas as pd
from snakemake.io import expand
from snakemake.workflow import srcdir


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


def read_samples(samples_file, data_dir, strict_mode=False, lineages_csv=False):
    samples_file = find_input_file(samples_file)
    if '.tsv' in samples_file or '.csv' in samples_file:
        separator = '\t'
        if '.csv' in samples_file:
            separator = ','
        try:
            # signame col doesn't need to exist: "If you pass extra name in this list, it will add another new column with that name with NaN values"
            if not lineages_csv:
                samples = pd.read_csv(samples_file, dtype=str, sep=separator, names = ["sample", "filename", "signame"])
            else:
                samples = pd.read_csv(samples_file, dtype=str, sep=separator, header=0)
        except Exception as e:
            sys.stderr.write(f"\n\tError: {samples_file} file is not properly formatted. Please fix.\n\n")
            print(e)
        # if signame column was not given, it will be NaNs. Fill NA's with sample names.
        samples['signame'] = samples['signame'].fillna(samples['sample'])
    elif '.xls' in samples_file:
        try:
            if not lineages_csv:
                samples = pd.read_excel(samples_file, dtype=str, names = ["sample", "filename", "signame"])
            else:
                samples = pd.read_excel(samples_file, dtype=str, header=0)
        except Exception as e:
            sys.stderr.write(f"\n\tError: {samples_file} file is not properly formatted. Please fix.\n\n")
            print(e)
        # if signame column was not given, it will be NaNs. Fill NA's with sample names.
        samples['signame'] = samples['signame'].fillna(samples['sample'])
    else:
        sample_list = [ line.strip() for line in open(samples_file, 'rt') ]
        sample_list = [ line for line in sample_list if line ]   # remove empty lines
        samplesD = {"sample":sample_list,"filename":sample_list,"signame":sample_list} # maybe later try removing *fa.gz or the like
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


def sanitize_path(path):
    # expand `~`, get absolute path
    path = os.path.expanduser(path)
    path = os.path.abspath(path)
    return path


def make_db_fullname(row):
    row["db_fullname"] = row["db_basename"] + "." + row["alphabet"] + "-k" +  str(row["ksize"]) + "-scaled" + str(row["scaled"])
    return row


def load_database_info(databases_file, existing_db_info=None,strict_mode=False):
    db_file = find_input_file(databases_file)
    if '.tsv' in db_file or '.csv' in db_file:
        separator = '\t'
        if '.csv' in db_file:
            separator = ','
        try:
            db_info = pd.read_csv(db_file, sep=separator, index_col=False)
            assert list(db_info.columns) == ["db_basename","alphabet","ksize","scaled","path","info_path"]
            db_info = db_info.apply(make_db_fullname, axis=1)
            # verify_integrity checks the new index for duplicates.
            db_info.set_index("db_fullname", inplace=True, verify_integrity=True)
            if existing_db_info:
                db_info = pd.concat(existing_db_info, db_info, verify_integrity=True)
        except Exception as e:
            sys.stderr.write(f"\n\tError: {db_file} file is not properly formatted. Please fix.\n\n")
            print(e)
    elif '.xls' in db_file:
        try:
            db_info = pd.read_excel(db_file, index_col=False)
            assert list(db_info.columns) == ["db_basename","alphabet","ksize","scaled","path","info_path"]
            db_info = db_info.apply(make_db_fullname, axis=1)
            db_info.set_index("db_fullname", inplace=True, verify_integrity=True)
            if existing_db_info:
                db_info = pd.concat(existing_db_info, db_info, verify_integrity=True)
        except Exception as e:
            sys.stderr.write(f"\n\tError: {db_file} file is not a properly formatted excel file. Please fix.\n\n")
            print(e)
    return db_info


def check_dbinfo(config, db_basename, dbinfo, alphabet_info, db_input="default", strict_mode=False):
    databases_to_use=[]
    db_targs, database_targets, info_targets=[],[],[]
    db_suffix = config["database_suffix"]
    info_suffix = config["database_info_suffix"]

    # get database fullnames
    available_databases = list(dbinfo.index)

    # construct fullnames that we want to make
    desired_dbfullnames = []
    for alpha, alphaInfo in alphabet_info.items():
        desired_dbfullnames += expand("{name}.{alpha}-k{ksize}-scaled{scaled}", name=db_basename, alpha=alpha, ksize=alphaInfo["ksizes"], scaled=alphaInfo["scaled"])

    # this now checks alphas, ksizes, scaled!
    for db in set(desired_dbfullnames):
        if db in available_databases:
            databases_to_use.append(db)
            # get database_targets
            database_targets.append(f"{db}.{db_suffix}")
            info_targets.append(f"{db_basename}.{info_suffix}")
        else:
            # get alpha, ksize, scaled from db_fullname
            alpha_k_scaled = db.split(f"{db_basename}.")[1]
            this_alpha, this_k, this_scaled = alpha_k_scaled.split("-")
            print(f'** ERROR: A {db_input} {db_basename} database is not provided for {this_alpha}, {this_k}, {this_scaled}.')
            if strict_mode:
                print('** Strict mode is on. Exiting.')
                sys.exit(-1)
            else:
                print(f'Strict mode is off: attempting to continue. Removing {db} from search databases list.')
    # remove any duplicates
    info_targets=list(set(info_targets))
    db_targs = database_targets + info_targets
    return databases_to_use, db_targs

def check_and_set_pipeline_databases(config, db_info, strict_mode=False):
    databases, database_targets=[],[]

    # find the default databases for this pipeline
    pipeline = config["pipeline"]
    default_databases = config["pipelines"][pipeline].get("databases", [])
    # use pipeline dbs unless user turns them off
    no_use_defaults = config.get("turn_off_default_databases", False)

    alphabet_info = config["alphabet_info"]
    if not no_use_defaults:
        # check that we have the info for these default databases
        for db in default_databases:
            db_names, db_targets = check_dbinfo(config, db, db_info, alphabet_info, strict_mode=strict_mode)
            databases+=db_names
            database_targets+=db_targets

    ## add user databases if info is available
    user_dbs = config.get("search_databases", [])

    for db in user_dbs:
        # add database to list of databases to search
        db_names, db_targets = check_dbinfo(config, db, db_info, alphabet_info, db_input="user", strict_mode=strict_mode)
        databases+=db_names
        database_targets+=db_targets

    if not databases:
        avail_dbs = list(db_info.index)
        print('\n** ERROR: no valid databases selected. ')
        print('Please choose from the available databases '\
              'or do not disable default databases. \n' )
        print('Available databases are: ',*avail_dbs, sep="\n" )
        print('** exiting.')
        sys.exit(-1)

    config["databases"] = databases

    # sanitize database directory path
    database_dir = sanitize_path(config["database_dir"])
    config["database_dir"] = database_dir
    final_db_targs = [os.path.join(database_dir, x) for x in database_targets]
    config["database_info"] = db_info

    return final_db_targs

def integrate_user_config(config):
    # first, check and set alphabets
    check_and_set_alphabets(config)
    # get info for pipeline we're running
    pipeline = config["pipeline"]
    db_info = ""
    database_targets=[]
    if config["pipelines"][pipeline]["databases_required"]:
       # if using databases, load database info and check compatibility with alphabets
        config["get_databases"] = True
        # get information for available databases (ksizes, file dl info, etc)
        default_dbinfo = config["default_database_info"]
        default_db_file = srcdir(default_dbinfo)
        db_info = load_database_info(default_db_file)
        # integrate user database info
        if config.get("user_database_info"):
            db_info = load_database_info(config["user_database_info"], existing_db_info=db_info)

        database_targets = check_and_set_pipeline_databases(config, db_info)
    else:
        config["get_databases"] = False

    return database_targets


def generate_database_targets(config):
    # if just downloading databases, allow folks to use minimal config
    if not config.get("input_type"):
        # if don't specify, use nucl to download all dbs
        config["input_type"] = "nucleotide"
    database_targets = integrate_user_config(config)
    return database_targets

def generate_targets(config, samples, output_dir="", generate_db_targets=False):
    pipeline_targets=[]
    ## integrate user settings and inputs
    database_targets = integrate_user_config(config)
    ## What alphabets are we using? ##
    alphabet_info = config["alphabet_info"]
    ## set run basename
    basename = config["basename"]

    ## What databases are we using?
    databases = config.get("databases", [])

    # TO DO: do this above with fullnames? Also make a rule to print csv info for generated indices
    index_names = []
    pipeline = config["pipeline"]
    #if pipeline == "generate_index" or pipeline == "compare_samples": eh. just always make index names, in case we need em
    for alpha, alphaInfo in alphabet_info.items():
        index_names+=expand("{basename}.{alpha}-k{ksize}-scaled{scaled}", basename=basename, alpha=alpha, ksize=alphaInfo["ksizes"], scaled=alphaInfo["scaled"])

    # generate targets for each step
    steps = config["pipelines"][pipeline]["steps"]
    for step in steps:
        step_outdir = config[step]["output_dir"]
        step_files = config[step]["output_files"]

        # fill variables in the output filenames
        for stepF in step_files:
            pipeline_targets += expand(os.path.join(output_dir, step_outdir, stepF), sample=samples, database=config.get("databases", []), basename=basename, index=index_names)

    if generate_db_targets:
        targets = database_targets + pipeline_targets
        return targets

    return pipeline_targets

