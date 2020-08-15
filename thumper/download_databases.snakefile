# download databases

import os, sys
import thumper.utils as tp

# should carry over from main snakefile.
# keeping here in case we want to use as standalone snake
database_dir = config['database_dir']
db_logs = os.path.join(database_dir, "logs")
db_benchmarks = os.path.join(database_dir, "benchmarks")

urls_begin = ["http", "ftp"]

database_info=config["databases"]


# if running as standalone, use this as rule all
rule download_databases:
    input: tp.generate_database_targets(config)

rule get_dbinfo:
    output:
        os.path.join(database_dir, "{db_name}-info.csv")
    params:
         csv_info= lambda w: database_info[w.db_name]["info_csv"]
    log: os.path.join(db_logs, "get_dbs", "{db_name}.info.get")
    run:
        if params.csv_info.startswith(tuple(urls_begin)):
            shell("curl -L {params.csv_info}  > {output}")
        else:
            full_input = os.path.abspath(str(params.csv_info))
            full_output = os.path.abspath(str(output))
            shell("ln -s {full_input} {full_output} 2> {log}")

rule get_sbt:
    output:
        os.path.join(database_dir, "{db_name}.{alphabet}-k{ksize}.sbt.zip")
    params:
        #sbt_info = lambda w: database_links[w.db_name]["sbt"]
        sbt_info = lambda w: database_info[w.db_name]["alphabets"][w.alphabet]["k" + w.ksize]["sbt"]
    log: os.path.join(db_logs, "get_dbs", "{db_name}.{alphabet}-k{ksize}.sbt.zip.get")
    run:
        if params.sbt_info.startswith(tuple(urls_begin)):
            shell("curl -L {params.sbt_info}  > {output}")
        else:
            full_input = os.path.abspath(str(params.sbt_info))
            full_output = os.path.abspath(str(output))
            shell("ln -s {full_input} {full_output} 2> {log}")

rule get_lca:
    output:
        os.path.join(database_dir, "{db_name}.{alphabet}.k{ksize}.lca.json.gz")
    params:
        lca_info = lambda w: database_info[w.db_name]["alphabets"][w.alphabet]["k" + w.ksize]["lca"]
    log: os.path.join(db_logs, "get_dbs", "{db_name}.{alphabet}-k{ksize}.lca.json.gz.get")
    run:
        if params.lca_info.startswith(tuple(urls_begin)):
            shell("curl -L {params.lca_info}  > {output}")
        else:
            full_input = os.path.abspath(str(params.lca_info))
            full_output = os.path.abspath(str(output))
            shell("ln -s {full_input} {full_output} 2> {log}")
