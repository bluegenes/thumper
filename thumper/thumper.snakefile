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
basename = config["basename"]

# check strict and force values
strict_val = config.get('strict', '1')
strict_mode = int(strict_val)
if not strict_mode:
    print('** WARNING: strict mode is OFF. Config errors will not force exit.')

force = config.get('force', '0')
force = int(force)
force_param = ''
if force:
    force_param = '--force'

## integrate values from user config ##
tp.check_and_set_alphabets(config, strict_mode=strict_mode)
alphabet_info = config["alphabet_info"]

if config.get("sample_list"):
    sample_info = tp.read_samples(config["sample_list"], data_dir)
else:
    print('** Error: Please provide proteomes/genomes as a txt file ' \
        '(one filename per line) or a csv file("sample,filename"; no headers ' \
        'using "sample_list:" in the config')
    sys.exit(-1)

sample_names = sample_info.index.tolist()

# snakemake info and workflow #
ascii_thumper = srcdir("animals/thumper")
failwhale = srcdir("animals/failwhale")

wildcard_constraints:
    alphabet="protein|dayhoff|hp|nucleotide", #|dna|rna",
    ksize="\d+",
    sample="\w+"
    #database = "(?!x\.).+"

onstart:
    print("------------------------------")
    print("Perform taxonomic classification using protein k-mers")
    print("------------------------------")


onsuccess:
    print("\n--- Workflow executed successfully! ---\n")
    shell('cat {ascii_thumper}')

onerror:
    print("Alas!\n")
    shell('cat {failwhale}')


# targeting rules
rule download_databases:
    input:
        expand(os.path.join(database_dir, "{database}.sbt.zip"), database=tp.check_databases(config))

rule mag_taxonomy:
    input: 
        expand(os.path.join(out_dir, "classify", "{basename}.x.{database}.taxonomy-report.csv"), basename=basename, database=tp.check_databases(config)),
        expand(os.path.join(out_dir, "classify", "{basename}.x.{database}.charcoal-lineages.csv"), basename=basename, database=tp.check_databases(config))
        #expand(os.path.join(out_dir, "classify", "{basename}.x.{database}.taxonomy-report.csv", basename=basename, database=config["databases"])
        #expand(os.path.join(out_dir, "classify", "{basename}.x.{database}.charcoal-lineages.csv"), basename=basename, database=config["databases"])
        
rule mag_contig_taxonomy:
    input:
        #expand(os.path.join(out_dir, "classify", "{basename}.x.{database}.taxonomy-report.csv", basename=basename, database=config["databases"])
        expand(os.path.join(out_dir, "classify", "{basename}.x.{database}.taxonomy-report.csv"), basename=basename, database=tp.check_databases(config))
    
#rule metagenome_classify:

rule build_index:
    input: 
        expand(os.path.join(out_dir, "index", "{index}.sbt.zip"), index = tp.build_index_names(config))
        #expand(os.path.join(out_dir, "index", "{index}.sbt.zip", index = config["index_names"])

include: "get_databases.snakefile"
include: "index.snakefile"
include: "common.snakefile"

## individual rules and steps ##

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
        gather_json=os.path.join(out_dir, 'contig-search', '{sample}.x.{db_name}.{alphabet}-k{ksize}-scaled{scaled}.contigs.gather.json'),
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
        gather_json=os.path.join(out_dir, 'genome-search', '{sample}.x.{db_name}.{alphabet}-k{ksize}-scaled{scaled}.gather.json'),
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

rule aggregate_genome_lca_gather:
    input: expand(os.path.join(out_dir, 'genome-search', '{sample}.x.{{db_name}}.{{alphabet}}-k{{ksize}}-scaled{{scaled}}.rankgather.csv'), sample=sample_names),
    output: os.path.join(out_dir, 'genome-search', f"{basename}" + ".x.{db_name}.{alphabet}-k{ksize}-scaled{scaled}.rankgather.csv")
    run:
        with open(str(output), 'w') as outF:
            write_header=True
            for inF in input:
                with open(str(inF), "r") as infile:
                    if write_header:
                        outF.write(infile.read())
                        write_header=False
                    else:
                        next(infile)
                        outF.write(infile.read())


# for a summary across all genomes, same database:
#def aggregate_taxonomy_files_by_database(w):
#    sample_list = config["sample_list"]
#    db_info = config["database_info"][w.db_name]["info_csv"],
#    json_tax,search_sigs=[],[]
#    tax_file = "classify/{sample}.x.{{db_name}}.{alphabet}-k{ksize}-scaled{scaled}.contigs-tax.json"
#    sig_file = "search/{sample}.x.{{db_name}}.{alphabet}-k{ksize}-scaled{scaled}.matches.sig"
#    for alpha, alphaInfo in alphabet_info.items():
#        json_tax += expand(os.path.join(out_dir, tax_file), sample=sample_names, alphabet=alpha, ksize = alphaInfo["ksizes"])
#        search_sigs += expand(os.path.join(out_dir, sig_file), sample=sample_names, alphabet=alpha, ksize = alphaInfo["ksizes"])
#
#    db_files = {"sample_list": sample_list,
#                "db_info": db_info,
#                "all_json": json_tax,
#                "all_sig": search_sigs}
#    return db_files

# papermill reporting rules
rule set_kernel:
    output:
        f"{out_dir}/.kernel.set"
    conda: 'envs/reporting-env.yml'
    resources:
        mem_mb=lambda wildcards, attempt: attempt *1000,
        runtime=1200,
    shell:
        """
        python -m ipykernel install --user --name thumper
        touch {output}
        """

localrules: make_genome_notebook, make_index, set_kernel, aggregate_gather_resultfiles

rule aggregate_gather_json:
    input:
        genome=expand(os.path.join(out_dir, 'genome-search', '{sample}.x.{{database}}.gather.json'), sample=sample_names),
    output:
        os.path.join(out_dir, "genome-search", "{basename}.x.{database}.gather.txt")
    resources:
        mem_mb=lambda wildcards, attempt: attempt *1000,
        runtime=200,
    run:
        with open(str(output), "w") as outF:
            header = ["name", "database", "genome-gather"]
            outF.write(",".join(header) + "\n")
            for sample in sample_names:
                genome_json = os.path.join('genome-search', f'{sample}.x.{wildcards.database}.gather.json')
                outF.write(sample + "," + wildcards.database + ',' + genome_json + ',' + "\n")

rule aggregate_contig_gather_json:
    input:
        contigs=expand(os.path.join(out_dir, 'contig-search', '{sample}.x.{{database}}.contigs.gather.json'), sample=sample_names),
    output:
        os.path.join(out_dir, "contig-search", "{basename}.x.{database}.gather.txt")
    resources:
        mem_mb=lambda wildcards, attempt: attempt *1000,
        runtime=200,
    run:
        with open(str(output), "w") as outF:
            header = ["name", "database", "contig-gather"]
            outF.write(",".join(header) + "\n")
            for sample in sample_names:
                contig_json = os.path.join('contig-search', f'{sample}.x.{wildcards.database}.contigs.gather.json')
                outF.write(sample + "," + wildcards.database + ',' + contig_json + "\n")

rule taxonomy_report:
    input:
        genome_info=os.path.join(out_dir, "genome-search", "{basename}.x.{database}.gather.txt"),
        #contig_info=os.path.join(out_dir, "contig-search", "{basename}.x.{database}.gather.txt"),
    output:
        genome_report=os.path.join(out_dir, "classify", "{basename}.x.{database}.taxonomy-report.csv"),
        charcoal_lineages=os.path.join(out_dir, "classify", "{basename}.x.{database}.charcoal-lineages.csv"),
        #common_contamination=os.path.join(out_dir, "classify", "{basename}.x.{database}.contamination-summary.json"),
        #contig_details_summary=os.path.join(out_dir, "classify", "{basename}.x.{database}.contig-details-summary.csv"),
    log: os.path.join(logs_dir, "genome-report", "{basename}.x.{database}.genome-report.log")
    benchmark: os.path.join(benchmarks_dir, "genome-report", "{basename}.x.{database}.genome-report.benchmark")
    resources:
        mem_mb=lambda wildcards, attempt: attempt *5000,
        runtime=6000,
    shell:
        """
        python -m thumper.compare_taxonomy \
                  --jsoninfo-file {input.gather_info} \
                  --match-rank "genus" \
                  --gather_min_matches 3 \
                  --min_f_ident 0.1 \
                  --min_f_major 0.2 \
                  --lineages-for-charcoal {output.charcoal_lineages} \
                  --output-csv {output.genome_report} 2> {log}
        """
        #          --contam-summary-json {output.common_contamination} \
        #          --contig-details-summary {output.contig_details_summary} 2> {log}
        #"""

rule report_genome_lineage:
    input:
        genome=expand(os.path.join(out_dir, 'genome-search', '{{sample}}.x.{database}.gather.json'), database=config['databases']),
    output:
        expand(os.path.join(report_dir, "{basename}.taxonomy-report.csv"), basename=basename)
    resources:
        mem_mb=lambda wildcards, attempt: attempt *1000,
        runtime=1200,
    shell:
        """
        python -m thumper.genome_report {input} -o {output}
        """


# use this notebook to aggregate files from 1. search containment, gather, 2. multiple alphas, 3. multiple databases?
rule make_genome_notebook:
    input:
        nb='thumper/notebooks/genome-report.ipynb',
        contigs=expand(os.path.join(out_dir, 'contig-search', '{{sample}}.x.{database}.contigs.gather.json'), database=config['databases']),
        genome=expand(os.path.join(out_dir, 'genome-search', '{{sample}}.x.{database}.rankgather.csv'), database=config['databases']),
        kernel_set = rules.set_kernel.output,
    params:
        name = lambda w: f"{w.sample}",
        databases= ",".join(config["databases"]),
        directory = os.path.abspath(out_dir),
    output:
        nb=os.path.join(report_dir, '{sample}.fig.ipynb'),
        html=os.path.join(report_dir, "{sample}.fig.html")
    conda: 'envs/reporting-env.yml'
    resources:
        mem_mb=lambda wildcards, attempt: attempt *1000,
        runtime=1200,
    shell:
        """
        papermill {input.nb} - -k thumper --cwd {report_dir} \
              -p directory {params.directory:q} -p render '' \
              -p name {params.name:q} \
              -p databases {params.databases:q} \
              > {output.nb}
        python -m nbconvert {output.nb} --to html --stdout --no-input --ExecutePreprocessor.kernel_name=thumper > {output.html}
        """

rule make_index:
    input:
        notebook='thumper/notebooks/report-index.ipynb',
        summary=os.path.join(out_dir, "classify", "{basename}.x.{database}.taxonomy-report.csv"),
        kernel_set = rules.set_kernel.output
    output:
        nb=os.path.join(report_dir, "{basename}.x.{database}.index.ipynb"),
        html=os.path.join(report_dir, "{basename}.x.{database}.index.html")
    params:
        directory = os.path.abspath(out_dir),
    conda: 'envs/reporting-env.yml'
    shell:
        """
        papermill {input.notebook} - -p name {wildcards.basename:q} -p render '' \
            -p database {wildcards.database:q} -p directory {params.directory:q} \
            -k thumper --cwd {report_dir} > {output.nb}
        python -m nbconvert {output.nb} --to html --stdout --no-input > {output.html}
        """
