## borrows heavily from charcoal, dammit, sgc cli infra
"Enable python -m thumper"
import os
import sys
import yaml
import glob
import subprocess

import click

def get_snakefile_path(name):
    thisdir = os.path.dirname(__file__)
    snakefile = os.path.join(thisdir, name)
    return snakefile


def get_package_configfile(filename):
    thisdir = os.path.dirname(__file__)
    configfile = os.path.join(thisdir, '.config', filename)
    return configfile


def run_snakemake(configfile, no_use_conda=False, no_use_mamba=False, verbose=False,
                  snakefile_name='classify.snakefile', extra_args=[]):
    # find the Snakefile relative to package path
    snakefile = get_snakefile_path(snakefile_name)

    # basic command
    cmd = ["snakemake", "-s", snakefile]

    # add --use-conda
    if not no_use_conda:
        if not no_use_mamba:
            cmd += ["--use-conda", "--conda-frontend", "mamba"]
        else:
            cmd += ["--use-conda"]

    # snakemake sometimes seems to want a default -j; set it to 1 for now.
    # can overridden later on command line.
    cmd += ["-j", "1"]

    # add rest of snakemake arguments
    cmd += list(extra_args)

    if configfile:
        # add defaults and system config files, in that order
        configfiles = [get_package_configfile("config.yaml"),
                       get_package_configfile("databases.yaml"),
                       get_package_configfile("pipelines.yaml"),
                       configfile]


        cmd += ["--configfile"] + configfiles

    if verbose:
        print('final command:', cmd)


    # runme
    try:
        subprocess.check_call(cmd)
    except subprocess.CalledProcessError as e:
        print(f'Error in snakemake invocation: {e}', file=sys.stderr)
        return e.returncode

#
# actual command line functions
#

@click.group()
def cli():
    pass

# create a run subcommand that by default passes all of its arguments
# on to snakemake (after setting Snakefile and config)
@click.command(context_settings={"ignore_unknown_options": True})
@click.argument('configfile')
@click.option('--no-use-conda', is_flag=True, default=False)
@click.option('--no-use-mamba', is_flag=True, default=False)
@click.option('--verbose', is_flag=True)
@click.argument('snakemake_args', nargs=-1)
def run(configfile, snakemake_args, no_use_conda, no_use_mamba, verbose):
    "execute thumper workflow (using snakemake underneath)"
    run_snakemake(configfile, snakefile_name='classify.snakefile',
                  no_use_conda=no_use_conda, no_use_mamba=no_use_mamba, verbose=verbose,
                  extra_args=snakemake_args)

# download databases using a special Snakefile
@click.command()
def download_db():
    "download the necessary databases"
    run_snakemake(None, snakefile_name='download_databases.snakefile',
                  no_use_conda=True)

# 'check' command
@click.command()
@click.argument('configfile')
def check(configfile):
    "check configuration"
    run_snakemake(configfile, extra_args=['check'])

# 'showconf' command
@click.command()
@click.argument('configfile')
def showconf(configfile):
    "show full configuration across default, system and project config files"
    run_snakemake(configfile, extra_args=['showconf'])

# 'info' command
@click.command()
def info():
    "provide basic install/config file info"
    from .version import version
    print(f"""
This is thumper version v{version}

Package install path: {os.path.dirname(__file__)}
Install-wide config file: {get_package_configfile('system.conf')}
snakemake Snakefile: {get_snakefile_path('Snakefile')}
""")

# 'init' command
@click.command()
@click.argument('configfile')
@click.option('--data-dir', nargs=1)
#@click.option('--lineages', nargs=1, default="")
@click.option('-f', '--force', is_flag=True)
def init(configfile, data_dir, lineages, force):
    "create a new, empty config file."
    stubname = os.path.basename(configfile)
    if configfile.endswith('.yaml'):
        stubname = stubname[:-5]
    else:
        configfile += '.yaml'

    if os.path.exists(configfile) and not force:
        print(f"** ERROR: configfile '{configfile}' already exists.")
        return -1

    sample_list = f'{stubname}.sample-list.txt'
    if data_dir:
        if os.path.exists(sample_list) and not force:
            print(f"** ERROR: sample list file '{sample_list}' already exists.")
            return -1
        samples = glob.glob(f'{data_dir}/*.fa')
        samples += glob.glob(f'{data_dir}/*.fna')
        samples += glob.glob(f'{data_dir}/*.faa')
        samples += glob.glob(f'{data_dir}/*.pep')
        samples += glob.glob(f'{data_dir}/*.fa.gz')
        samples += glob.glob(f'{data_dir}/*.fna.gz')
        samples += glob.glob(f'{data_dir}/*.faa.gz')
        samples += glob.glob(f'{data_dir}/*.pep.gz')
        print(f'found {len(samples)} samples in {data_dir}/*.{{fa,fna,faa,pep}}{{,.gz}}')
        samples = [ os.path.basename(g) for g in samples ]
        with open(sample_list, 'wt') as fp:
            fp.write("\n".join(samples))
        print(f"created '{sample_list}' with {len(samples)} samples in it.")

    #if lineages:
    #    print(f"Using provided lineages from '{lineages}'")
    #else:
    #    print("(No provided lineages file given.)")

    print(f"creating configfile '{configfile}' for project '{stubname}'")
    with open(configfile, 'wt') as fp:
        fp.write(\
f"""\
# location for all generated files
output_dir: output.{stubname}/

# list of sample filenames to classify
sample_list: {stubname}.sample-list.txt

# directory in which sample filenames live
data_dir: {data_dir}

# (optional) list of lineages for input samples. comment out or leave
# blank if none.
provided_lineages: {lineages}

# match_rank is the rank _above_ which cross-lineage matches are considered
# contamination. e.g. if set to 'superkingdom', then Archaeal matches in
# Bacterial samples will be contamination, but nothing else.
#
# values can be superkingdom, phylum, class, order, family, or genus.
match_rank: order
""")

cli.add_command(run)
cli.add_command(check)
cli.add_command(showconf)
cli.add_command(info)
cli.add_command(init)
cli.add_command(download_db)

def main():
    cli()

if __name__ == '__main__':
    main()