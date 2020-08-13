import os, sys
import numpy as np
import pandas as pd
import glob

from snakemake.workflow import srcdir

out_dir = config["thumper_dir"]
logs_dir = os.path.join(out_dir, "logs")
benchmarks_dir = os.path.join(out_dir, "benchmarks")
database_dir = config['database_dir']

#wildcard_constraints:
    #size="\d+"
    #transcriptome = config['transcriptome_name'],
    #database = "(?!x\.).+"

# snakemake workflow
onstart:
    print("------------------------------")
    print("Perform taxonomic classification using protein k-mers")
    print("------------------------------")

#thumper = srcdir("utils/animals/thumper")
onsuccess:
    print("\n--- Workflow executed successfully! ---\n")
    shell('cat {srcdir("utils/animals/thumper")}')

onerror:
    print("Alas!\n")

# include the databases snake
include: "download_databases.snakefile"

def make_sketch_type(input_type):
    if input_type == 'nucleotide':
        args = "translate"
    else:
        args = "protein"

    ## get the desired ksize, scaled values from the config
    #ksize = config["alphabet"][input_type]
    return args

#def make_sketch_args(signature_params):
#    
#    return args


rule sourmash_sketch:
    input: os.path.join(data_dir, "{sample}.sig")
    output: os.path.join(output_dir, "{sketchtype}", "{sample}.sig") # IF doing both nucleotide, protein --> need alphabet wildcard in here
    params:
        k= lambda w: (int(w.k) * ksize_multiplier[w.alphabet]),
        sketch_cmd = make_sketch_args(config['input_type']) 
        signame = lambda w: accession2signame[w.accession],
        abund_cmd = "--track-abundance",
    threads: 1
    resources:
        mem_mb=lambda wildcards, attempt: attempt *1000,
        runtime=1200,
    log: os.path.join(logs_dir, "sourmash", "{accession}_{alphabet}_scaled{scaled}_k{k}.dna.compute.log")
    benchmark: os.path.join(benchmarks_dir, "sourmash", "{accession}_{alphabet}_scaled{scaled}_k{k}.dna.compute.benchmark")
    conda: "envs/sourmash-dev.yml"
    shell:
        """
        sourmash sketch -k {params.k} --scaled={wildcards.scaled}  \
        {input} -o {output} {params.alpha_cmd} {params.abund_cmd} --merge={params.signame:q} 2> {log}
        """

# gather each sig

## enable multiple databases ---> impacts naming (cant just do x.__database__.)
rule sourmash_gather:
    input:
        query= rules.sourmash_sketch.output
        databases=config["gather_database"] 
    output:
        csv = os.path.join(out_dir, "gather", "{sample}.gather-matches.csv"),
        matches = os.path.join(out_dir, "gather", "{sample}.gather-matches.sig"),
        txt = os.path.join(out_dir, "gather", "{sample}.gather-matches.txt"),
        unassigned = os.path.join(out_dir, "gather", "{sample}.gather-unassigned.sig"),
    params:
        alpha_cmd = lambda w: "--" + config["alphabet"], # this just enables one ... if we want nucl + prot, change here
        scaled = config["gather_scaled"],
    resources:
        mem_mb=lambda wildcards, attempt: attempt *10000,
        runtime=600,
    log: os.path.join(logs_dir, "gather", "{sample}.gather.log")
    benchmark: os.path.join(benchmarks_dir, "gather", "{sample}.gather.benchmark")
    conda: "envs/sourmash-dev.yml"
    shell:
        # --ignore-abundance to turn abund off
        """
        sourmash gather {input.query} {input.databases} -o {output.csv} \
        {params.alpha_cmd} --scaled {params.scaled} \
        --save-matches {output.matches} --threshold-bp=0  \
        --output-unassigned {output.unassigned} \
        >& {output.txt} 2> {log}
        touch {output.csv} {output.matches}
        """
# touch empty output file to enable rna ones to fail (to do: handle failures properly downstream)


rule gather_to_tax:
    input:
        #gather_csv = rules.gather_sig.output.csv,
        gather_csv = os.path.join(gather_dir, "{db_name}.{alphabet}-k{ksize}","{genome}_x_{db_name}.{alphabet}-k{ksize}.gather.csv"),
        #lineages_csv = lambda w: refInfo[w.db_name]["lineages_csv"]
        lineages_csv = lambda w: refInfo[w.alphabet][w.db_name]["lineages_csv"]
        #lineages_csv="/home/ntpierce/2020-gtdb-smash/gtdb-lineages.rna-filenames.n0th-representative-at-genus.csv"
    output:
        gather_tax = os.path.join(gather_dir, "{db_name}.{alphabet}-k{ksize}", "{genome}_x_{db_name}.{alphabet}-k{ksize}.gather_summary.csv"),
        top_matches = os.path.join(gather_dir, "{db_name}.{alphabet}-k{ksize}", "{genome}_x_{db_name}.{alphabet}-k{ksize}.gather_tophits.csv")
    log: os.path.join(logs_dir, "gather_to_tax", "{db_name}.{alphabet}-k{ksize}", input_type, "{genome}_x_{db_name}.{alphabet}-k{ksize}.gather-to-tax.log")
    benchmark: os.path.join(benchmarks_dir, "gather_to_tax", "{db_name}.{alphabet}-k{ksize}","{genome}_x_{db_name}.{alphabet}-k{ksize}.gather-to-tax.benchmark")
    group: "gather"
    resources:
        mem_mb=lambda wildcards, attempt: attempt *10000,
        runtime=200,
    conda: "envs/sourmash-dev.yml"
    shell:
        """
        python scripts/gather-to-tax.py {input.gather_csv} {input.lineages_csv} --tophits-csv {output.top_matches} > {output.gather_tax} 2> {log}
        """


##### traverse directory is troublesome here --> each set of gather requires its own folder
rule aggregate_gather_to_tax:
    # make spreadsheet: each proteome:: top lineage hit
    input:
        gather_tophits= lambda w: expand(os.path.join(gather_dir, "{{db_name}}.{{alphabet}}-k{{ksize}}", "{genome}_x_{{db_name}}.{{alphabet}}-k{{ksize}}.gather_tophits.csv"), genome=  sampleInfo[w.sample]["accessions"])
    output:
        summary_csv=os.path.join(summary_dir, "{sample}_x_{db_name}.{alphabet}-k{ksize}.gather_tophits.csv"),
    params:
        gather_dir= os.path.join(gather_dir, "{db_name}.{alphabet}-k{ksize}"),
        true_lineages_cmd = lambda w: "--true-lineages-csv " + sampleInfo[w.sample].get("true_lineages_csv") if sampleInfo[w.sample].get("true_lineages_csv") else "",
    #log: os.path.join(logs_dir, "aggregate_gather", "{sample}_x_{db_name}." + f"{alphabet}-k{ksize}.gather_tophits.log")
    log: os.path.join(logs_dir, "aggregate_gather", "{sample}_x_{db_name}.{alphabet}-k{ksize}.gather_tophits.log")
    benchmark: os.path.join(benchmarks_dir, "aggregate_gather", "{sample}_x_{db_name}.{alphabet}-k{ksize}.gather_tophits.benchmark")
    resources:
        mem_mb=lambda wildcards, attempt: attempt *10000,
        runtime=60,
    conda: "envs/sourmash-dev.yml"
    shell:
        #python scripts/aggregate-gather-to-tax-tophits.py --input-is-directory --output-csv {output} {params.gather_dir} --true-lineages-csv {input.true_lineages} 2> {log}
        """
        python scripts/aggregate-gather-to-tax-tophits.py --input-is-directory --output-csv {output} {params.gather_dir} {params.true_lineages_cmd} 2> {log}
        """
