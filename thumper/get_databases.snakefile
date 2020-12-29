# download databases

import os, sys
import thumper.utils as tp

# should carry over from main snakefile.
# keeping here in case we want to use as standalone snake
database_dir = config['database_dir']
db_logs = os.path.join(database_dir, "logs")
db_benchmarks = os.path.join(database_dir, "benchmarks")

urls_begin = ["http", "ftp"]

database_info=config["database_info"]


localrules: get_dbinfo, get_sbt

rule get_dbinfo:
    output:
        os.path.join(database_dir, "{db_basename}.info.csv")
    params:
         #csv_info= lambda w: database_info.at[w.db_name, 'info_path']
         csv_info= lambda w: config["database_info"].loc[database_info["db_basename"]== w.db_basename]["info_path"][0]
    log: os.path.join(db_logs, "get_dbs", "{db_basename}.info.get")
    threads: 1
    resources:
        mem_mb=1000,
        runtime=600
    run:
        if params.csv_info.startswith(tuple(urls_begin)):
            shell("curl -L {params.csv_info}  > {output}")
        else:
            full_input = os.path.abspath(str(params.csv_info))
            full_output = os.path.abspath(str(output))
            shell("ln -s {full_input} {full_output} 2> {log}")

rule get_sbt:
    output:
        os.path.join(database_dir, "{database}.sbt.zip")
    params:
        #sbt_info = lambda w: database_info[w.db_name]["alphabets"][w.alphabet]["k" + str(w.ksize)]["sbt"]
        sbt_info= lambda w: config["database_info"].at[w.database,"path"]
    log: os.path.join(db_logs, "get_dbs", "{database}.sbt.get")
    threads: 1
    resources:
        mem_mb=1000,
        runtime=600
    run:
        if params.sbt_info.startswith(tuple(urls_begin)):
            shell("curl -L {params.sbt_info}  > {output}")
        else:
            full_input = os.path.abspath(str(params.sbt_info))
            full_output = os.path.abspath(str(output))
            shell("ln -s {full_input} {full_output} 2> {log}")


# if running as standalone, use this as rule all
#rule download_databases:
#    input: tp.generate_database_targets(config)

#rule get_lca:
#    output:
#        os.path.join(database_dir, "{db_name}.{alphabet}.k{ksize}.lca.json.gz")
#    params:
#        lca_info = lambda w: database_info[w.db_name]["alphabets"][w.alphabet]["k" + str(w.ksize)]["lca"]
#    log: os.path.join(db_logs, "get_dbs", "{db_name}.{alphabet}-k{ksize}.lca.json.gz.get")
#    threads: 1
#    resources:
#        mem_mb=1000,
#        runtime=600
#    run:
#        if params.lca_info.startswith(tuple(urls_begin)):
#            shell("curl -L {params.lca_info}  > {output}")
#        else:
#            full_input = os.path.abspath(str(params.lca_info))
#            full_output = os.path.abspath(str(output))
#            shell("ln -s {full_input} {full_output} 2> {log}")
