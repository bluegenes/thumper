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
report_dir = os.path.join(out_dir, "reports")


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

ascii_thumper = srcdir("animals/thumper")
failwhale = srcdir("animals/failwhale")
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
    include: "get_databases.snakefile"
include: "index.snakefile"
include: "common.snakefile"
alphabet_info = config["alphabet_info"]

def build_sketch_params(output_type):
    sketch_cmd = ""
    input_type = config["input_type"]
    # if input is dna, build dna, translate sketches
    if input_type == "nucleotide":
        if output_type == "nucleotide":
            ksizes = config["alphabet_info"]["nucleotide"].get("ksizes", config["alphabet_defaults"]["nucleotide"]["ksizes"])
            scaled = config["alphabet_info"]["nucleotide"].get("scaled", config["alphabet_defaults"]["nucleotide"]["scaled"])
            # always track abund when sketching (?)
            sketch_cmd = "dna -p " + "k=" + ",k=".join(map(str, ksizes)) + f",scaled={str(scaled)}" + ",abund"
            return sketch_cmd
        else:
            sketch_cmd = "translate "
    else:
        # if input is protein, just build protein sketches
        sketch_cmd = "protein "
    for alpha in ["protein", "dayhoff", "hp"]:
        if alpha in config["alphabet_info"].keys():
        ## default build protein, dayhoff, hp sigs at the default ksizes from config
            ksizes = config["alphabet_info"][alpha].get("ksizes", config["alphabet_defaults"][alpha]["ksizes"])
            scaled = config["alphabet_info"][alpha].get("scaled", config["alphabet_defaults"][alpha]["scaled"])
        else:
            ksizes = config["alphabet_defaults"][alpha]["ksizes"]
            scaled = config["alphabet_defaults"][alpha]["scaled"]
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
        database=os.path.join(database_dir, "{db_name}.{alphabet}-k{ksize}-scaled{scaled}.sbt.zip")
    output:
        csv = os.path.join(out_dir, "search", "{sample}.x.{db_name}.{alphabet}-k{ksize}-scaled{scaled}.matches.csv"),
        matches = os.path.join(out_dir, "search", "{sample}.x.{db_name}.{alphabet}-k{ksize}-scaled{scaled}.matches.sig"),
        txt = os.path.join(out_dir, "search", "{sample}.x.{db_name}.{alphabet}-k{ksize}-scaled{scaled}.matches.txt"),
    params:
        alpha_cmd = lambda w: alphabet_info[w.alphabet]["alpha_cmd"],
        scaled = lambda w: alphabet_info[w.alphabet]["scaled"],
        ksize = lambda w: int(w.ksize)*int(alphabet_info[w.alphabet]["ksize_multiplier"]),
        search_threshold = float(config.get("search_threshold", 0.001))
    resources:
        mem_mb=lambda wildcards, attempt: attempt *20000,
        runtime=6000,
    log: os.path.join(logs_dir, "search", "{sample}.x.{db_name}.{alphabet}-k{ksize}-scaled{scaled}.search.log")
    benchmark: os.path.join(benchmarks_dir, "search", "{sample}.x.{db_name}.{alphabet}-k{ksize}-scaled{scaled}.search.benchmark")
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
rule contigs_search:
    input:
        sample_file=lambda w: os.path.join(data_dir, sample_info.at[w.sample, 'filename']),
        sample_sig=os.path.join(out_dir, "signatures", "{sample}.sig"),
        matches=os.path.join(out_dir, "search", "{sample}.x.{db_name}.{alphabet}-k{ksize}-scaled{scaled}.matches.sig"),
        db_info=lambda w: os.path.join(database_dir, f"{w.db_name}.info.csv")
    output:
        search_csv=os.path.join(out_dir, 'contig-search', '{sample}.x.{db_name}.{alphabet}-k{ksize}-scaled{scaled}.contigs.search.csv'),
        search_sig=os.path.join(out_dir, 'contig-search', '{sample}.x.{db_name}.{alphabet}-k{ksize}-scaled{scaled}.contigs.search.matches.sig'),
        ranksearch_csv=os.path.join(out_dir, 'contig-search', '{sample}.x.{db_name}.{alphabet}-k{ksize}-scaled{scaled}.contigs.ranksearch.csv'),
        ranksearch_sig=os.path.join(out_dir, 'contig-search', '{sample}.x.{db_name}.{alphabet}-k{ksize}-scaled{scaled}.contigs.ranksearch.matches.sig'),
        rankgather_csv=os.path.join(out_dir, 'contig-search', '{sample}.x.{db_name}.{alphabet}-k{ksize}-scaled{scaled}.contigs.rankgather.csv'),
        unmatched=os.path.join(out_dir, 'contig-search', '{sample}.x.{db_name}.{alphabet}-k{ksize}-scaled{scaled}.contigs.unmatched.fq'),
    params:
        ksize = lambda w: int(w.ksize)*int(alphabet_info[w.alphabet]["ksize_multiplier"]),
        moltype = lambda w: alphabet_info[w.alphabet]["moltype"],
        out_prefix = lambda w: os.path.join(out_dir, 'contig-search', f"{w.sample}.x.{w.db_name}.{w.alphabet}-k{w.ksize}-scaled{w.scaled}"),
    conda: 'envs/sourmash-dev.yml'
    resources:
        mem_mb=lambda wildcards, attempt: attempt *100000,
        runtime=6000,
    log: os.path.join(logs_dir, "contig-search", "{sample}.x.{db_name}.{alphabet}-k{ksize}-scaled{scaled}.contigs-search.log")
    benchmark: os.path.join(benchmarks_dir, "contig-search", "{sample}.x.{db_name}.{alphabet}-k{ksize}-scaled{scaled}.contigs-search.benchmark")
    shell: 
        """
        python -m thumper.search_or_gather \
            --genome {input.sample_file} \
            --genome-sig {input.sample_sig} \
            --matches-sig {input.matches} \
            --lineages-csv {input.db_info} \
            --alphabet {params.moltype} \
            --ksize {params.ksize} \
            --gather \
            --output-prefix {params.out_prefix}
        """

rule genome_search:
    input:
        sample_file=lambda w: os.path.join(data_dir, sample_info.at[w.sample, 'filename']),
        sample_sig=os.path.join(out_dir, "signatures", "{sample}.sig"),
        matches=os.path.join(out_dir, "search", "{sample}.x.{db_name}.{alphabet}-k{ksize}-scaled{scaled}.matches.sig"),
        #db_info=config["database_info"].loc["gtdb-nine.hp-k42-scaled10"]["info_path"]
        db_info=lambda w: os.path.join(database_dir, f"{w.db_name}.info.csv")
    output:
        search_csv=os.path.join(out_dir, 'genome-search', '{sample}.x.{db_name}.{alphabet}-k{ksize}-scaled{scaled}.search.csv'),
        search_sig=os.path.join(out_dir, 'genome-search', '{sample}.x.{db_name}.{alphabet}-k{ksize}-scaled{scaled}.search.matches.sig'),
        ranksearch_csv=os.path.join(out_dir, 'genome-search', '{sample}.x.{db_name}.{alphabet}-k{ksize}-scaled{scaled}.ranksearch.csv'),
        ranksearch_sig=os.path.join(out_dir, 'genome-search', '{sample}.x.{db_name}.{alphabet}-k{ksize}-scaled{scaled}.ranksearch.matches.sig'),
        rankgather_csv=os.path.join(out_dir, 'genome-search', '{sample}.x.{db_name}.{alphabet}-k{ksize}-scaled{scaled}.rankgather.csv'),
    params:
        ksize = lambda w: int(w.ksize)*int(alphabet_info[w.alphabet]["ksize_multiplier"]),
        moltype = lambda w: alphabet_info[w.alphabet]["moltype"],
        out_prefix = lambda w: os.path.join(out_dir, 'genome-search', f"{w.sample}.x.{w.db_name}.{w.alphabet}-k{w.ksize}-scaled{w.scaled}"),
    conda: 'envs/sourmash-dev.yml'
    resources:
        mem_mb=lambda wildcards, attempt: attempt *100000,
        runtime=6000,
    log: os.path.join(logs_dir, "genome-search", "{sample}.x.{db_name}.{alphabet}-k{ksize}-scaled{scaled}.genome-search.log")
    benchmark: os.path.join(benchmarks_dir, "genome-search", "{sample}.x.{db_name}.{alphabet}-k{ksize}-scaled{scaled}.genome-search.benchmark")
    shell:
        """
        python -m thumper.search_or_gather \
            --genome {input.sample_file} \
            --genome-sig {input.sample_sig} \
            --matches-sig {input.matches} \
            --lineages-csv {input.db_info} \
            --alphabet {params.moltype} \
            --ksize {params.ksize} \
            --search-genome --gather --no-search-contigs \
            --output-prefix {params.out_prefix}
        """

# for a summary across all genomes, same database:
def aggregate_taxonomy_files_by_database(w):
    sample_list = config["sample_list"]
    db_info = config["database_info"][w.db_name]["info_csv"],
    json_tax,search_sigs=[],[]
    tax_file = "classify/{sample}.x.{{db_name}}.{alphabet}-k{ksize}-scaled{scaled}.contigs-tax.json"
    sig_file = "search/{sample}.x.{{db_name}}.{alphabet}-k{ksize}-scaled{scaled}.matches.sig"
    for alpha, alphaInfo in alphabet_info.items():
        json_tax += expand(os.path.join(out_dir, tax_file), sample=sample_names, alphabet=alpha, ksize = alphaInfo["ksizes"])
        search_sigs += expand(os.path.join(out_dir, sig_file), sample=sample_names, alphabet=alpha, ksize = alphaInfo["ksizes"])

    db_files = {"sample_list": sample_list,
                "db_info": db_info,
                "all_json": json_tax,
                "all_sig": search_sigs}
    return db_files


# papermill reporting rules
rule set_kernel:
    output:
        f"{out_dir}/.kernel.set"
    conda: 'envs/reporting-env.yml'
    shell: """
        python -m ipykernel install --user --name thumper
        touch {output}
    """

# use this notebook to aggregate files from 1. search containment, gather, 2. multiple alphas, 3. multiple databases?
rule make_notebook:
    input:
        nb='thumper/notebooks/genome-report.ipynb',
        contigs=os.path.join(out_dir, 'genome-search', '{sample}.x.{db_name}.{alphabet}-k{ksize}-scaled{scaled}.contigs.rankgather.csv'),
        genome=os.path.join(out_dir, 'genome-search', '{sample}.x.{db_name}.{alphabet}-k{ksize}-scaled{scaled}.rankgather.csv'),
        kernel_set = rules.set_kernel.output,
    params:
        name = lambda w: f"{w.sample}.x.{w.db_name}.{w.alphabet}-k{w.ksize}-scaled{w.scaled}"
    output:
        os.path.join(report_dir, 'genome-search', '{sample}.x.{db_name}.{alphabet}-k{ksize}-scaled{scaled}.fig.ipynb')
    conda: 'envs/reporting-env.yml'
    shell: 
        """
        papermill {input.nb} - -k thumper --cwd {report_dir} \
              -p directory .. -p render '' \
              -p name {params.name:q} \
              > {output}
        """

#rule make_html:
#    input:
#        notebook=os.path.join(report_dir, 'genome-search', '{sample}.x.{db_name}.{alphabet}-k{ksize}-scaled{scaled}.fig.ipynb')
#        genome_summary=f'{out_dir}/genome_summary.csv',
#        hitlist=f'{out_dir}/hit_list_for_filtering.csv',
#        contigs_json=f'{out_dir}/{{g}}.contigs-tax.json',
#    output:
#        report_dir + '/{g}.fig.html',
#    conda: 'envs/reporting-env.yml'
#    shell: 
#        """
#        python -m nbconvert {input.notebook} --stdout --no-input --ExecutePreprocessor.kernel_name=thumper > {output}
#        """
     

#rule make_index:
#    input:
#        notebook='thumper/notebooks/report_index.ipynb',
#        summary=f'{output_dir}/genome_summary.csv',
#        kernel_set = rules.set_kernel.output
#    output:
#        nb=f'{report_dir}/index.ipynb',
#        html=f'{report_dir}/index.html',
#    conda: 'envs/reporting-env.yml'
#    shell: 
#        """
#        papermill {input.notebook} - -p name {out_dir:q} -p render '' \
#            -p directory .. -k thumper --cwd {report_dir} > {output.nb}
#        python -m nbconvert {output.nb} --stdout --no-input > {output.html}
#        """
