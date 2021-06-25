import os
import sys
import pandas as pd
from snakemake.io import expand
from snakemake.workflow import srcdir

__all__ = ['find_input_file', 'read_samples',
           'check_and_set_alpha_ksize',
           'load_databases_csv', 'find_valid_databases',
           'check_and_set_databases', 'process_user_config']


def sanitize_path(path):
    # expand `~`, get absolute path
    path = os.path.expanduser(path)
    path = os.path.abspath(path)
    return path


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
    if not found_file:
        raise ValueError(f'Cannot find specified {name} file {filename}')
    if verbose:
        sys.stderr.write(f'\tFound {name} at {found_file}\n')
    return found_file


def _file_exists(filename, data_dir):
    fname = os.path.join(data_dir, filename)
    return os.path.isfile(fname)


def read_samples(config, *, strict=False):
    try:
        filename = config.get('sample_info', '')
        samples_file = find_input_file(filename)
    except ValueError as exc:
        print('** ERROR: Please provide proteomes/genomes as a txt file ' \
        '(one filename per line) or a csv file("sample,filename"; no headers ' \
        'using "sample_list:" in the config')
        sys.exit(-1)
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

    data_dir = config.get('data_dir', '')
    if data_dir:
        data_dir = sanitize_path(data_dir)
    # check if file exists
    samples['exists'] = samples['filename'].apply(_file_exists, data_dir=data_dir)
    if not pd.Series(samples['exists']).all():
        missing_files = samples[samples['exists'] == False].index
        if strict:
            raise ValueError(f'Sample files {",".join(missing_files)} do not exist in {data_dir}')
        else:
            print(f'Sample files {",".join(missing_files)} do not exist in {data_dir}')
            print('Strict mode is OFF. Attempting to continue.')
            samples = samples[samples['exists'] == True]
    samples.set_index("sample", inplace=True)
    config['samples'] = samples
    return config


def check_and_set_alpha_ksize(config, *, strict=False):
    """
    Check the input type, alphabet, and ksize specified in the configfile.
    Set alphabet and ksize for sketch and search.
    """
    alphabet_defaults = config["alphabet_defaults"]
    input_type = config.get("input_type")
    alphabet = config.get('alphabet')
    ksize = config.get('ksize')
    if not input_type:
        raise ValueError('"input_type: protein" or "input_type: nucleotide" must be specified in the config file')
    # nucleotide input
    elif input_type in ["nucleotide", "dna", "rna"]:
        if alphabet in ["nucleotide", "dna", "rna"]:
            config["sketch_type"] = 'nucleotide'
            if not ksize:
                config["ksize"] = alphabet_defaults["nucleotide"]["ksizes"]
            config["scaled"] = alphabet_defaults["nucleotide"]["scaled"]
        # translate to protein
        elif alphabet in ["protein", "dayhoff", "hp"]:
            config["sketch_type"] = protein
            if not ksize:
                config["ksize"] = alphabet_defaults["nucleotide"]["ksizes"]
            config["scaled"] = alphabet_defaults["nucleotide"]["scaled"]
    #protein input
    elif input_type in ["protein"]:
        config["sketch_type"] = "protein"
        if not ksize:
            config["ksize"] = alphabet_defaults[alphabet]["ksizes"]
        config["scaled"] = alphabet_defaults[alphabet]["scaled"]

    else: #unknown input
        raise ValueError(f'input type {input_type} must be "protein" or "nucleotide"')

    #check ksizes again; make sure is list
    ksizes = config['ksize']
    if not isinstance(ksize, list):
        config['ksize'] = [ksizes]
    return config


def _make_db_fullname(row):
    row["db_fullname"] = row["db_basename"] + "." + row["alphabet"] + "-k" +  str(row["ksize"])
    return row


def load_databases_csv(databases_file, existing_db_info=None):
    db_file = find_input_file(databases_file)
    db_info=None
    if '.tsv' in db_file or '.csv' in db_file:
        separator = '\t'
        if '.csv' in db_file:
            separator = ','
        try:
            db_info = pd.read_csv(db_file, sep=separator, index_col=False)
            assert list(db_info.columns) == ["db_basename","alphabet","ksize","scaled","path","taxonomy_path"]
            db_info = db_info.apply(_make_db_fullname, axis=1)
            # verify_integrity checks the new index for duplicates.
            db_info.set_index("db_fullname", inplace=True, verify_integrity=True)
            if existing_db_info:
                db_info = pd.concat(existing_db_info, db_info, verify_integrity=True)
        except Exception as e:
            print(e)
            raise ValueError(f"\n\tError: {db_file} file is not properly formatted. Please fix.\n\n")
    elif '.xls' in db_file:
        try:
            db_info = pd.read_excel(db_file, index_col=False)
            assert list(db_info.columns) == ["db_basename","alphabet","ksize","scaled","path","info_path"]
            db_info = db_info.apply(_make_db_fullname, axis=1)
            db_info.set_index("db_fullname", inplace=True, verify_integrity=True)
            if existing_db_info:
                db_info = pd.concat(existing_db_info, db_info, verify_integrity=True)
        except Exception as e:
            print(e)
            raise ValueError(f"\n\tError: {db_file} file is not a properly formatted excel file. Please fix.\n\n")
    return db_info


def find_valid_databases(databases, db_info, config, *, strict=False):
    databases_to_use=[]
    db_basenames = []
    alpha = config['alphabet']
    ksizes = config['ksize']
    for db in databases:
        for k in ksizes:
            #check that this alpha/ksize exists in our database info file
            db_fullname = f"{db}.{alpha}-k{k}"
            if db_fullname in db_info.index:
                db_basenames.append(db)
                databases_to_use.append(db_fullname)
            else:
                if strict:
                    raise ValueError(f'The {db} database is not provided for alphabet:{alpha}, ksize:{ksize}')
                else:
                    print(f'Strict mode is off: attempting to continue. Removing {db} from search databases list.')

    if not databases_to_use:
        raise ValueError('No valid databases provided.')

    return databases_to_use, db_basenames


def check_and_set_databases(config, *, strict=False):
    # first, get load all the database info from default and user CSVs
    default_db_file = srcdir(config["default_database_info"])
    db_info = load_databases_csv(default_db_file)
    # integrate user database info
    if config.get("user_database_info"):
       # load user database info
       db_info = load_databases_csv(config["user_database_info"], existing_db_info=db_info)
    # store database info in config to enable easy download and use
    config["database_info"] = db_info

    #now check if specified search databases are valid
    databases = config.get('search_databases', [])
    if not isinstance(databases, list):
        databases = [databases]
    try:
        valid_databases, valid_db_basenames= find_valid_databases(databases, db_info, config)
    except ValueError as exc:
        avail_dbs = list(db_info.index)
        print('Available databases are: ',*avail_dbs, sep="\n" )
        raise
    config['valid_databases'] = valid_databases
    config['valid_db_basenames'] = valid_db_basenames
    return config


def process_user_config(config):
    strict_mode = config.get('strict_mode', False)
    try:
        samples = read_samples(config, strict=strict_mode)
        config = check_and_set_alpha_ksize(config, strict=strict_mode)
        config = check_and_set_databases(config, strict=strict_mode)
    except ValueError as exc:
        print(f'ERROR: Could not properly integrate configuration values.')
        print(f'{str(exc)}')
        sys.exit(-1)

    return config

