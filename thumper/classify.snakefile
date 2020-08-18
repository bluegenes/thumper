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

wildcard_constraints:
    prot_alphabet="protein|dayhoff|hp",
    nucl_alphabet="nucleotide|dna|rna",
    ksize="\d+"
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
    input: tp.generate_targets(config, sample_list, out_dir, generate_db_targets=False)

# include the databases, common utility snakefiles
include: "download_databases.snakefile"
include: "common.snakefile"

def build_sketch_params(output_type):
    # todo: name alphabet_defaults better to enable user config override (or addition??) 
    sketch_cmd = ""
    input_type = config["input_type"]
    # if input is dna, build dna, translate sketches
    if input_type == "nucleotide":
        if output_type == "nucleotide":
            ksizes = config["alphabet_defaults"]["nucleotide"]["ksizes"]
            scaled = config["alphabet_defaults"]["nucleotide"]["scaled"]
            # always track abund when sketching
            sketch_cmd = "-p " + "k=" + ",k=".join(map(str, ksizes)) + f",scaled={str(scaled)}" + ",abund"
            #print(sketch_cmd)
            return sketch_cmd
        else:
            sketch_cmd = "translate "
    else:
        # if input is protein, just build protein sketches
        sketch_cmd = "protein "
    for alpha in ["protein", "dayhoff", "hp"]:
        ## default build protein, dayhoff, hp sigs at the default ksizes from config
        ksizes = config["alphabet_defaults"][alpha]["ksizes"]
        scaled = config["alphabet_defaults"][alpha]["scaled"]
        sketch_cmd += " -p " + alpha + ",k=" + ",k=".join(map(str, ksizes)) + f",scaled={str(scaled)}" + ",abund"
    #print(sketch_cmd)
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

## if want to enable multiple databases in same cmd, will need to change naming /targeting scheme
rule sourmash_gather_protein:
    input:
        prot_query=rules.sourmash_sketch_protein.output,
        database=os.path.join(database_dir, "{db_name}.{prot_alphabet}-k{ksize}.sbt.zip")
    output:
        csv = os.path.join(out_dir, "gather", "{filename}.x.{db_name}.{prot_alphabet}-k{ksize}.gather-matches.csv"),
        matches = os.path.join(out_dir, "gather", "{filename}.x.{db_name}.{prot_alphabet}-k{ksize}.gather-matches.sig"),
        txt = os.path.join(out_dir, "gather", "{filename}.x.{db_name}.{prot_alphabet}-k{ksize}.gather-matches.txt"),
        unassigned = os.path.join(out_dir, "gather", "{filename}.x.{db_name}.{prot_alphabet}-k{ksize}.gather-unassigned.sig"),
    params:
        alpha_cmd = "", #"--protein", #lambda w: "--" + config["alphabet"], # one alpha at at time. rn using ksize instead
        scaled = config["gather_scaled"],
        ksize = lambda w: int(w.ksize)*3
    resources:
        mem_mb=lambda wildcards, attempt: attempt *10000,
        runtime=6000,
    log: os.path.join(logs_dir, "gather", "{filename}.x.{db_name}.{prot_alphabet}-k{ksize}.gather.log")
    benchmark: os.path.join(benchmarks_dir, "gather", "{filename}.x.{db_name}.{prot_alphabet}-k{ksize}.gather.benchmark")
    conda: "envs/sourmash-dev.yml"
    shell:
        # --ignore-abundance to turn abund off --> enable this in params?
        """
        sourmash gather {input.prot_query} {input.database} -o {output.csv} \
        {params.alpha_cmd} --scaled {params.scaled} \
        -k {params.ksize}  \
        --save-matches {output.matches} --threshold-bp=0  \
        --output-unassigned {output.unassigned} \
        >& {output.txt} 2> {log}
        touch {output.csv} {output.matches} {output.unassigned}
        """

rule sourmash_search_containment_protein:
    input:
        prot_query=rules.sourmash_sketch_protein.output,
        database=os.path.join(database_dir, "{db_name}.{prot_alphabet}-k{ksize}.sbt.zip")
    output:
        csv = os.path.join(out_dir, "search-containment", "{filename}.x.{db_name}.{prot_alphabet}-k{ksize}.search-contain-matches.csv"),
        matches = os.path.join(out_dir, "search-containment", "{filename}.x.{db_name}.{prot_alphabet}-k{ksize}.search-contain-matches.sig"),
        txt = os.path.join(out_dir, "search-containment", "{filename}.x.{db_name}.{prot_alphabet}-k{ksize}.search-contain-matches.txt"),
    params:
        alpha_cmd = "", #"--protein", #lambda w: "--" + config["alphabet"], # one alpha at at time. rn using ksize instead
        scaled = config["gather_scaled"],
        ksize = lambda w: int(w.ksize)*3
    resources:
        mem_mb=lambda wildcards, attempt: attempt *10000,
        runtime=6000,
    log: os.path.join(logs_dir, "search-contain", "{filename}.x.{db_name}.{prot_alphabet}-k{ksize}.search-contain.log")
    benchmark: os.path.join(benchmarks_dir, "search-contain", "{filename}.x.{db_name}.{prot_alphabet}-k{ksize}.search-contain.benchmark")
    conda: "envs/sourmash-dev.yml"
    shell:
        # --ignore-abundance to turn abund off --> enable this in params?
        """
        sourmash search --containment {input.prot_query} {input.database} \
        -o {output.csv} \
        {params.alpha_cmd} --scaled {params.scaled} \
        -k {params.ksize}  --threshold=0.001 \
        --save-matches {output.matches}  \
        >& {output.txt} 2> {log}
        touch {output.csv} {output.matches}
        """


rule sourmash_gather_nucleotide:
    input:
        nucl_query=rules.sourmash_sketch_nucleotide.output,
        database=os.path.join(database_dir, "{db_name}.{nucl_alphabet}-k{ksize}.sbt.zip")
    output:
        csv = os.path.join(out_dir, "gather", "{filename}.x.{db_name}.{nucl_alphabet}-k{ksize}.gather-matches.csv"),
        matches = os.path.join(out_dir, "gather", "{filename}.x.{db_name}.{nucl_alphabet}-k{ksize}.gather-matches.sig"),
        txt = os.path.join(out_dir, "gather", "{filename}.x.{db_name}.{nucl_alphabet}-k{ksize}.gather-matches.txt"),
        unassigned = os.path.join(out_dir, "gather", "{filename}.x.{db_name}.{nucl_alphabet}-k{ksize}.gather-unassigned.sig"),
    params:
        alpha_cmd = "", #"--protein", #lambda w: "--" + config["alphabet"], # one alpha at at time
        scaled = config["gather_scaled"],
    resources:
        mem_mb=lambda wildcards, attempt: attempt *10000,
        runtime=6000,
    log: os.path.join(logs_dir, "gather", "{filename}.x.{db_name}.{nucl_alphabet}-k{ksize}.gather.log")
    benchmark: os.path.join(benchmarks_dir, "gather", "{filename}.x.{db_name}.{nucl_alphabet}-k{ksize}.gather.benchmark")
    conda: "envs/sourmash-dev.yml"
    shell:
        # --ignore-abundance to turn abund off --> enable this in params?
        """
        sourmash gather {input.nucl_query} {input.database} -o {output.csv} \
        {params.alpha_cmd} --scaled {params.scaled} \
        -k {wildcards.ksize}  \
        --save-matches {output.matches} --threshold-bp=0  \
        --output-unassigned {output.unassigned} \
        >& {output.txt} 2> {log}
        touch {output.csv} {output.matches}
        """


rule sourmash_search_containment_nucleotide:
    input:
        nucl_query=rules.sourmash_sketch_nucleotide.output,
        database=os.path.join(database_dir, "{db_name}.{nucl_alphabet}-k{ksize}.sbt.zip")
    output:
        csv = os.path.join(out_dir, "search-containment", "{filename}.x.{db_name}.{nucl_alphabet}-k{ksize}.search-contain-matches.csv"),
        matches = os.path.join(out_dir, "search-containment", "{filename}.x.{db_name}.{nucl_alphabet}-k{ksize}.search-contain-matches.sig"),
        txt = os.path.join(out_dir, "search-containment", "{filename}.x.{db_name}.{nucl_alphabet}-k{ksize}.search-contain-matches.txt"),
    params:
        alpha_cmd = "", #"--protein", #lambda w: "--" + config["alphabet"], # one alpha at at time. rn using ksize instead
        scaled = config["gather_scaled"],
    resources:
        mem_mb=lambda wildcards, attempt: attempt *10000,
        runtime=6000,
    log: os.path.join(logs_dir, "search-contain", "{filename}.x.{db_name}.{nucl_alphabet}-k{ksize}.search-contain.log")
    benchmark: os.path.join(benchmarks_dir, "search-contain", "{filename}.x.{db_name}.{nucl_alphabet}-k{ksize}.search-contain.benchmark")
    conda: "envs/sourmash-dev.yml"
    shell:
        # --ignore-abundance to turn abund off --> enable this in params?
        """
        sourmash search --containment {input.nucl_query} {input.database} \
        -o {output.csv} --threshold=0.001 \
        {params.alpha_cmd} --scaled {params.scaled} \
        -k {wildcards.ksize}  \
        --save-matches {output.matches}  \
        >& {output.txt} 2> {log}
        touch {output.csv} {output.matches}
        """

#rule contigs_clean_just_taxonomy:
rule contig_classify:
    input: 
        sample_file=os.path.join(data_dir, "{filename}"),
        #matches=rules.sourmash_search_containment_protein.output.matches,
        matches=os.path.join(out_dir, "search-containment", "{filename}.x.{db_name}.{prot_alphabet}-k{ksize}.search-contain-matches.sig"),
        db_info=lambda w: config["databases"][w.db_name]["info_csv"],
#        script = srcdir('just_taxonomy.py'),
#        genome = genome_dir + '/{f}',
#        matches = output_dir + '/{f}.gather-matches.sig',
#        lineages = config['lineages_csv']
    output: 
        clean=os.path.join(out_dir,"classify", "{filename}.x.{db_name}.{prot_alphabet}-k{ksize}.clean.fa.gz"),
        dirty=os.path.join(out_dir,"classify", "{filename}.x.{db_name}.{prot_alphabet}-k{ksize}.dirty.fa.gz"),
        report=os.path.join(out_dir,"classify", "{filename}.x.{db_name}.{prot_alphabet}-k{ksize}.report.txt"),
        contig_report=os.path.join(out_dir, "classify", "{filename}.x.{db_name}.{prot_alphabet}-k{ksize}.contigs.csv"),
        csv=os.path.join(out_dir, "classify", "{filename}.x.{db_name}.{prot_alphabet}-k{ksize}.summary.csv")
    params:
        lineage = "", #get_provided_lineage,
        force = "", #force_param,
        match_rank = "genus", #match_rank,
        #moltype = config['moltype']
        moltype = lambda w: w.prot_alphabet,
    conda: "envs/sourmash-dev.yml"
    shell: """
        python -m thumper.charcoal_just_taxonomy \
            --genome {input.sample_file} --lineages_csv {input.db_info} \
            --matches_sig {input.matches} \
            --clean {output.clean} --dirty {output.dirty} \
            --report {output.report} --summary {output.csv} \
            --match-rank {params.match_rank} \
            --contig-report {output.contig_report}
    """





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
