import os, sys
import pandas as pd
import glob

from snakemake.workflow import srcdir
import thumper.utils as tp

out_dir = config["output_dir"]
logs_dir = os.path.join(out_dir, "logs")
benchmarks_dir = os.path.join(out_dir, "benchmarks")
database_dir = config['database_dir']
data_dir = config['data_dir'].rstrip('/')


strict_val = config.get('strict', '1')
strict_mode = int(strict_val)
if not strict_mode:
    print('** WARNING: strict mode is OFF. Config errors will not force exit.')

force = config.get('force', '0')
force = int(force)
force_param = ''
if force:
    force_param = '--force'


# snakemake workflow

#wildcard_constraints:
    #size="\d+"
    #database = "(?!x\.).+"

sample_list = tp.read_samples(config["sample_list"], data_dir)

onstart:
    print("------------------------------")
    print("Perform taxonomic classification using protein k-mers")
    print("------------------------------")

ascii_thumper = srcdir("utils/animals/thumper")
failwhale = srcdir("utils/animals/failwhale")
onsuccess:
    print("\n--- Workflow executed successfully! ---\n")
    shell('cat {ascii_thumper}')

onerror:
    print("Alas!\n")
    shell('cat {failwhale}')

rule all:
    input: tp.generate_targets(config, sample_list, out_dir)

# include the databases, common utility snakefiles
include: "download_databases.snakefile"
include: "common.snakefile"

## TODO - sketch params are built differently!
# k=31,scaled=1000,abund
#-p k=20,num=500,protein -p k=19,num=400,dayhoff,abund
def build_sketch_params(output_type):
    ### if input is dna, build dna, translate sketches
    #### if input is protein, just build protein sketches
    ## default build protein, dayhoff, hp sigs at the default ksizes from config
    # always track abund when sketching
    ## turn this off with config
    # todo: name alphabet_defaults better to enable user config override (or addition??) 
    sketch_cmd = ""
    # build sketch cmd
    input_type = config["input_type"]
    if input_type == "nucleotide":
        if output_type == "nucleotide":
            ksizes = config["alphabet_defaults"]["nucleotide"]["ksizes"]
            scaled = config["alphabet_defaults"]["nucleotide"]["scaled"]
            sketch_cmd = "-p " + "k=" + ",k=".join(map(str, ksizes)) + f",scaled={str(scaled)}" + ",abund"
            print(sketch_cmd)
            return sketch_cmd
        else:
            sketch_cmd = "translate "
    else:
        sketch_cmd = "protein "
    for alpha in ["protein", "dayhoff", "hp"]:
        ksizes = config["alphabet_defaults"][alpha]["ksizes"]
        scaled = config["alphabet_defaults"][alpha]["scaled"]
        sketch_cmd += " -p " + alpha + ",k=" + ",k=".join(map(str, ksizes)) + f",scaled={str(scaled)}" + ",abund"
    print(sketch_cmd)
    return sketch_cmd



rule sourmash_sketch_nucleotide:
    input: os.path.join(data_dir, "{filename}")
    output:
        os.path.join(out_dir, "signatures", "{filename}.nucleotide.sig"),
    params:
        sketch_params = build_sketch_params("nucleotide"),
        #signame = lambda w: accession2signame[w.accession],
    threads: 1
    resources:
        mem_mb=lambda wildcards, attempt: attempt *1000,
        runtime=1200,
    log: os.path.join(logs_dir, "sourmash_sketch_nucl", "{filename}.nucl.log")
    benchmark: os.path.join(benchmarks_dir, "sourmash_sketch_nucl", "{filename}.nucl.benchmark")
    conda: "envs/sourmash-dev.yml"
    shell:
        """
        sourmash sketch dna {params.sketch_params} -o {output} {input}  2> {log}
        """

rule sourmash_sketch_protein:
    input: os.path.join(data_dir, "{filename}")
    output:
        os.path.join(out_dir, "signatures", "{filename}.protein.sig"),
    params:
        sketch_params = build_sketch_params("protein")
        #signame = lambda w: accession2signame[w.accession],
    threads: 1
    resources:
        mem_mb=lambda wildcards, attempt: attempt *1000,
        runtime=1200,
    log: os.path.join(logs_dir, "sourmash_sketch_prot", "{filename}.protein.log")
    benchmark: os.path.join(benchmarks_dir, "sourmash_sketch_prot", "{filename}.protein.benchmark")
    conda: "envs/sourmash-dev.yml"
    shell:
        """
        sourmash sketch {params.sketch_params} {input} -o {output}  2> {log}
        """

# gather each sig

## enable multiple databases ---> impacts naming (cant just do x.__database__.)
#rule sourmash_gather:
#    input:
#        query= rules.sourmash_sketch.output
#        databases=config["gather_database"] 
#    output:
#        csv = os.path.join(out_dir, "gather", "{sample}.gather-matches.csv"),
#        matches = os.path.join(out_dir, "gather", "{sample}.gather-matches.sig"),
#        txt = os.path.join(out_dir, "gather", "{sample}.gather-matches.txt"),
#        unassigned = os.path.join(out_dir, "gather", "{sample}.gather-unassigned.sig"),
#    params:
#        alpha_cmd = lambda w: "--" + config["alphabet"], # this just enables one ... if we want nucl + prot, change here
#        scaled = config["gather_scaled"],
#    resources:
#        mem_mb=lambda wildcards, attempt: attempt *10000,
#        runtime=600,
#    log: os.path.join(logs_dir, "gather", "{sample}.gather.log")
#    benchmark: os.path.join(benchmarks_dir, "gather", "{sample}.gather.benchmark")
#    conda: "envs/sourmash-dev.yml"
#    shell:
#        # --ignore-abundance to turn abund off
#        """
#        sourmash gather {input.query} {input.databases} -o {output.csv} \
#        {params.alpha_cmd} --scaled {params.scaled} \
#        --save-matches {output.matches} --threshold-bp=0  \
#        --output-unassigned {output.unassigned} \
#        >& {output.txt} 2> {log}
##        touch {output.csv} {output.matches}
#        """
## touch empty output file to enable rna ones to fail (to do: handle failures properly downstream)


#rule gather_to_tax:
#    input:
        #gather_csv = rules.gather_sig.output.csv,
#        gather_csv = os.path.join(gather_dir, "{db_name}.{alphabet}-k{ksize}","{genome}_x_{db_name}.{alphabet}-k{ksize}.gather.csv"),
#        #lineages_csv = lambda w: refInfo[w.db_name]["lineages_csv"]
#        lineages_csv = lambda w: refInfo[w.alphabet][w.db_name]["lineages_csv"]
#        #lineages_csv="/home/ntpierce/2020-gtdb-smash/gtdb-lineages.rna-filenames.n0th-representative-at-genus.csv"
#    output:
#        gather_tax = os.path.join(gather_dir, "{db_name}.{alphabet}-k{ksize}", "{genome}_x_{db_name}.{alphabet}-k{ksize}.gather_summary.csv"),
#        top_matches = os.path.join(gather_dir, "{db_name}.{alphabet}-k{ksize}", "{genome}_x_{db_name}.{alphabet}-k{ksize}.gather_tophits.csv")
#    log: os.path.join(logs_dir, "gather_to_tax", "{db_name}.{alphabet}-k{ksize}", input_type, "{genome}_x_{db_name}.{alphabet}-k{ksize}.gather-to-tax.log")
#    benchmark: os.path.join(benchmarks_dir, "gather_to_tax", "{db_name}.{alphabet}-k{ksize}","{genome}_x_{db_name}.{alphabet}-k{ksize}.gather-to-tax.benchmark")
#    group: "gather"
#    resources:
#        mem_mb=lambda wildcards, attempt: attempt *10000,
#        runtime=200,
#    conda: "envs/sourmash-dev.yml"
#    shell:
#        """
#        python scripts/gather-to-tax.py {input.gather_csv} {input.lineages_csv} --tophits-csv {output.top_matches} > {output.gather_tax} 2> {log}
#        """


##### traverse directory is troublesome here --> each set of gather requires its own folder
#rule aggregate_gather_to_tax:
#    # make spreadsheet: each proteome:: top lineage hit
#    input:
#        gather_tophits= lambda w: expand(os.path.join(gather_dir, "{{db_name}}.{{alphabet}}-k{{ksize}}", "{genome}_x_{{db_name}}.{{alphabet}}-k{{ksize}}.gather_tophits.csv"), genome=  sampleInfo[w.sample]["accessions"])
#    output:
#        summary_csv=os.path.join(summary_dir, "{sample}_x_{db_name}.{alphabet}-k{ksize}.gather_tophits.csv"),
    #params:
#        gather_dir= os.path.join(gather_dir, "{db_name}.{alphabet}-k{ksize}"),
#        true_lineages_cmd = lambda w: "--true-lineages-csv " + sampleInfo[w.sample].get("true_lineages_csv") if sampleInfo[w.sample].get("true_lineages_csv") else "",
#    #log: os.path.join(logs_dir, "aggregate_gather", "{sample}_x_{db_name}." + f"{alphabet}-k{ksize}.gather_tophits.log")
#    log: os.path.join(logs_dir, "aggregate_gather", "{sample}_x_{db_name}.{alphabet}-k{ksize}.gather_tophits.log")
#    benchmark: os.path.join(benchmarks_dir, "aggregate_gather", "{sample}_x_{db_name}.{alphabet}-k{ksize}.gather_tophits.benchmark")
#    resources:
#        mem_mb=lambda wildcards, attempt: attempt *10000,
#        runtime=60,
#    conda: "envs/sourmash-dev.yml"
#    shell:
#        #python scripts/aggregate-gather-to-tax-tophits.py --input-is-directory --output-csv {output} {params.gather_dir} --true-lineages-csv {input.true_lineages} 2> {log}
#        """
#        python scripts/aggregate-gather-to-tax-tophits.py --input-is-directory --output-csv {output} {params.gather_dir} {params.true_lineages_cmd} 2> {log}
#        """
