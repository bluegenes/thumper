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
    alphabet="protein|dayhoff|hp|nucleotide", #|dna|rna",
    ksize="\d+"
    #database = "(?!x\.).+"

if config.get("sample_list"):
    sample_info = tp.read_samples(config["sample_list"], data_dir)
else:
    print('** Error: Please provide proteomes/genomes as a txt file ' \
        '(one filename per line) or a csv file("sample,filename"; no headers ' \
        'using "sample_list:" in the config')
    sys.exit(-1)

sample_names = sample_info.index.tolist()

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
    input: tp.generate_targets(config, sample_names, out_dir, generate_db_targets=False)

# include the databases, index, common utility snakefiles
if config["get_databases"]:
    include: "download_databases.snakefile"
include: "index.snakefile"
include: "common.snakefile"
alphabet_info = config["alphabet_info"]

def build_sketch_params(output_type):
    sketch_cmd = ""
    input_type = config["input_type"]
    # if input is dna, build dna, translate sketches
    if input_type == "nucleotide":
        if output_type == "nucleotide":
            ksizes = config["alphabet_info"]["nucleotide"]["ksizes"]
            scaled = config["alphabet_info"]["nucleotide"]["scaled"]
            # always track abund when sketching (?)
            sketch_cmd = "dna -p " + "k=" + ",k=".join(map(str, ksizes)) + f",scaled={str(scaled)}" + ",abund"
            return sketch_cmd
        else:
            sketch_cmd = "translate "
    else:
        # if input is protein, just build protein sketches
        sketch_cmd = "protein "
    for alpha in ["protein", "dayhoff", "hp"]:
        ## default build protein, dayhoff, hp sigs at the default ksizes from config
        ksizes = config["alphabet_info"][alpha]["ksizes"]
        scaled = config["alphabet_info"][alpha]["scaled"]
        sketch_cmd += " -p " + alpha + ",k=" + ",k=".join(map(str, ksizes)) + f",scaled={str(scaled)}" + ",abund"
    return sketch_cmd

if config["input_type"] == "nucleotide":
    rule sourmash_sketch_nucleotide_input:
        input: lambda w: os.path.join(data_dir, sample_info.at[w.sample, 'filename'])
        output:
            full_sketch=os.path.join(out_dir, "signatures", "{sample}.sig"),
        params:
            nucl_sketch_params = build_sketch_params("nucleotide"),
            translate_sketch_params = build_sketch_params("protein"),
            nucl_sketch=os.path.join(out_dir, "signatures", "{sample}.nucleotide.sig"),
            prot_sketch=os.path.join(out_dir, "signatures", "{sample}.translate.sig"),
            #signame = lambda w: accession2signame[w.accession],
        threads: 1
        resources:
            mem_mb=lambda wildcards, attempt: attempt *1000,
            runtime=1200,
        log: os.path.join(logs_dir, "sourmash_sketch_nucl_input", "{sample}.sketch.log")
        benchmark: os.path.join(benchmarks_dir, "sourmash_sketch_nucl_input", "{sample}.sketch.benchmark")
        conda: "envs/sourmash-dev.yml"
        shell:
            """
            sourmash sketch {params.nucl_sketch_params} -o {params.nucl_sketch} --name {wildcards.sample} {input}  2> {log}
            sourmash sketch {params.translate_sketch_params} -o {params.prot_sketch} --name {wildcards.sample} {input}  2>> {log}
            sourmash sig cat {params.nucl_sketch} {params.prot_sketch} -o {output.full_sketch} 2>> {log}
            rm {params.nucl_sketch}
            rm {params.prot_sketch}
            """
else:
    rule sourmash_sketch_protein_input:
        input: lambda w: os.path.join(data_dir, sample_info.at[w.sample, 'filename'])
        output:
            os.path.join(out_dir, "signatures", "{sample}.sig"),
        params:
            sketch_params = build_sketch_params("protein")
            #signame = lambda w: accession2signame[w.accession],
        threads: 1
        resources:
            mem_mb=lambda wildcards, attempt: attempt *1000,
            runtime=1200,
        log: os.path. join(logs_dir, "sourmash_sketch_prot_input", "{sample}.sketch.log")
        benchmark: os.path.join(benchmarks_dir, "sourmash_sketch_prot_input", "{sample}.sketch.benchmark")
        #wildcard_constraints:
        #    alphabet="protein|dayhoff|hp",
        conda: "envs/sourmash-dev.yml"
        shell:
            """
            sourmash sketch {params.sketch_params} -o {output} --name {wildcards.sample} {input} 2> {log}
            """

rule sourmash_search_containment:
    input:
        query=os.path.join(out_dir, "signatures", "{sample}.sig"),
        database=os.path.join(database_dir, "{db_name}.{alphabet}-k{ksize}.sbt.zip")
    output:
        csv = os.path.join(out_dir, "search", "{sample}.x.{db_name}.{alphabet}-k{ksize}.matches.csv"),
        matches = os.path.join(out_dir, "search", "{sample}.x.{db_name}.{alphabet}-k{ksize}.matches.sig"),
        txt = os.path.join(out_dir, "search", "{sample}.x.{db_name}.{alphabet}-k{ksize}.matches.txt"),
    params:
        alpha_cmd = lambda w: alphabet_info[w.alphabet]["alpha_cmd"],
        scaled = lambda w: alphabet_info[w.alphabet]["scaled"],
        ksize = lambda w: int(w.ksize)*int(alphabet_info[w.alphabet]["ksize_multiplier"]),
        search_threshold = config.get("search_threshold", 0.001)
    resources:
        mem_mb=lambda wildcards, attempt: attempt *20000,
        runtime=6000,
    log: os.path.join(logs_dir, "search", "{sample}.x.{db_name}.{alphabet}-k{ksize}.search.log")
    benchmark: os.path.join(benchmarks_dir, "search", "{sample}.x.{db_name}.{alphabet}-k{ksize}.search.benchmark")
    conda: "envs/sourmash-dev.yml"
    shell:
        """
        sourmash search --containment {input.query} \
        {input.database} -o {output.csv} \
        {params.alpha_cmd} --scaled {params.scaled} \
        -k {params.ksize}  --threshold={params.search_threshold} \
        --save-matches {output.matches}  \
        >& {output.txt} 2> {log}
        touch {output.csv} {output.matches}
        """

# generate contigs taxonomy
rule contigs_taxonomy:
    input:
        sample_file=lambda w: os.path.join(data_dir, sample_info.at[w.sample, 'filename']),
        sig=os.path.join(out_dir, "signatures", "{sample}.sig"),
        matches=os.path.join(out_dir, "search", "{sample}.x.{db_name}.{alphabet}-k{ksize}.matches.sig"),
        db_info=lambda w: config["database_info"][w.db_name]["info_csv"],
    output:
        json=os.path.join(out_dir, 'classify', '{sample}.x.{db_name}.{alphabet}-k{ksize}.contigs-tax.json'),
    params:
        ksize = lambda w: int(w.ksize)*int(alphabet_info[w.alphabet]["ksize_multiplier"]),
    conda: 'envs/sourmash-dev.yml'
    resources:
        mem_mb=lambda wildcards, attempt: attempt *100000,
        runtime=6000,
    log: os.path.join(logs_dir, "classify", "{sample}.x.{db_name}.{alphabet}-k{ksize}.contigs-tax.log")
    benchmark: os.path.join(benchmarks_dir, "classify", "{sample}.x.{db_name}.{alphabet}-k{ksize}.contigs-tax.benchmark")
    shell: 
        """
        python -m thumper.contigs_search \
            --genome {input.sample_file} \
            --lineages-csv {input.db_info} \
            --alphabet {wildcards.alphabet}\
            --ksize {params.ksize} \
            --genome-sig {input.sig} \
            --matches-sig {input.matches} \
            --json-out {output.json}
        """

# For later? Expand sample, alpha, ksize, do NOT expand dbname. Aggregate reports over alpha/ksize
def aggregate_taxonomy_files_by_database(w):
    sample_list = config["sample_list"]
    db_info = config["database_info"][w.db_name]["info_csv"],
    json_tax,search_sigs=[],[]
    tax_file = "classify/{sample}.x.{{db_name}}.{alphabet}-k{ksize}.contigs-tax.json"
    sig_file = "search/{sample}.x.{{db_name}}.{alphabet}-k{ksize}.matches.sig"
    for alpha, alphaInfo in alphabet_info.items():
        json_tax += expand(os.path.join(out_dir, tax_file), sample=sample_names, alphabet=alpha, ksize = alphaInfo["ksizes"])
        search_sigs += expand(os.path.join(out_dir, sig_file), sample=sample_names, alphabet=alpha, ksize = alphaInfo["ksizes"])

    db_files = {"sample_list": sample_list,
                "db_info": db_info,
                "all_json": json_tax,
                "all_sig": search_sigs}
    return db_files

rule make_hit_list:
    #input: unpack(aggregate_taxonomy_files)
    input:
        sample_list = config["sample_list"],
        db_info = lambda w: config["database_info"][w.db_name]["info_csv"],
        all_json = expand(os.path.join(out_dir, "classify/{sample}.x.{{db_name}}.{{alphabet}}-k{{ksize}}.contigs-tax.json"), sample=sample_names),
        all_sig = expand(os.path.join(out_dir, "search/{sample}.x.{{db_name}}.{{alphabet}}-k{{ksize}}.matches.sig"), sample=sample_names),
    output:
        os.path.join(out_dir, "classify", "{basename}.x.{db_name}.{alphabet}-k{ksize}.hit_list_for_filtering.csv")
    params:
        output_dir = out_dir,
        min_f_major = float(config["min_f_major"]),
        min_f_ident = float(config["min_f_ident"]),
        moltype = lambda w: alphabet_info[w.alphabet]["moltype"],
        ksize = lambda w: int(w.ksize)*int(alphabet_info[w.alphabet]["ksize_multiplier"]),
    resources:
        mem_mb=lambda wildcards, attempt: attempt *40000,
        runtime=6000,
    conda: 'envs/sourmash-dev.yml'
    shell:
        #    --provided-lineages {input.provided_lineages} \
        """
        python -m thumper.compare_taxonomy \
            --input-directory {params.output_dir} \
            --genome-list-file {input.sample_list} \
            --database-name {wildcards.db_name} \
            --lineages-csv {input.db_info} \
            --alphabet {params.moltype} \
            --ksize {params.ksize} \
            --output {output} \
            --min_f_ident={params.min_f_ident} \
            --min_f_major={params.min_f_major}
        """



