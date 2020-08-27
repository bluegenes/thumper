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

# include the databases, common utility snakefiles
if config["get_databases"]:
    include: "download_databases.snakefile"
include: "common.snakefile"

def build_sketch_params(output_type):
    sketch_cmd = ""
    input_type = config["input_type"]
    # if input is dna, build dna, translate sketches
    if input_type == "nucleotide":
        if output_type == "nucleotide":
            ksizes = config["alphabet_info"]["nucleotide"]["ksizes"]
            scaled = config["alphabet_info"]["nucleotide"]["scaled"]
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
        ksizes = config["alphabet_info"][alpha]["ksizes"]
        scaled = config["alphabet_info"][alpha]["scaled"]
        sketch_cmd += " -p " + alpha + ",k=" + ",k=".join(map(str, ksizes)) + f",scaled={str(scaled)}" + ",abund"
    #print(sketch_cmd)
    return sketch_cmd


rule sourmash_sketch_nucleotide:
    #input: os.path.join(data_dir, "{filename}")
    input: lambda w: os.path.join(data_dir, sample_info.at[w.sample, 'filename'])
    output:
        os.path.join(out_dir, "signatures", "{sample}.nucleotide.sig"),
    params:
        sketch_params = build_sketch_params("nucleotide"),
        #signame = lambda w: accession2signame[w.accession],
    threads: 1
    resources:
        mem_mb=lambda wildcards, attempt: attempt *1000,
        runtime=1200,
    log: os.path.join(logs_dir, "sourmash_sketch_nucl", "{sample}.nucl.log")
    benchmark: os.path.join(benchmarks_dir, "sourmash_sketch_nucl", "{sample}.nucl.benchmark")
    conda: "envs/sourmash-dev.yml"
    shell:
        """
        sourmash sketch dna {params.sketch_params} -o {output} --name {wildcards.sample} {input}  2> {log}
        """

rule sourmash_sketch_protein:
    input: lambda w: os.path.join(data_dir, sample_info.at[w.sample, 'filename'])
    #input: os.path.join(data_dir, "{filename}")
    output:
        os.path.join(out_dir, "signatures", "{sample}.protein.sig"),
    params:
        sketch_params = build_sketch_params("protein")
        #signame = lambda w: accession2signame[w.accession],
    threads: 1
    resources:
        mem_mb=lambda wildcards, attempt: attempt *1000,
        runtime=1200,
    log: os.path.join(logs_dir, "sourmash_sketch_prot", "{sample}.protein.log")
    benchmark: os.path.join(benchmarks_dir, "sourmash_sketch_prot", "{sample}.protein.benchmark")
    conda: "envs/sourmash-dev.yml"
    shell:
        """
        sourmash sketch {params.sketch_params} -o {output} --name {wildcards.sample} {input} 2> {log}
        """

rule sourmash_search_containment_protein:
    input:
        prot_query=rules.sourmash_sketch_protein.output,
        database=os.path.join(database_dir, "{db_name}.{prot_alphabet}-k{ksize}.sbt.zip")
    output:
        csv = os.path.join(out_dir, "search-containment", "{sample}.x.{db_name}.{prot_alphabet}-k{ksize}.search-contain-matches.csv"),
        matches = os.path.join(out_dir, "search-containment", "{sample}.x.{db_name}.{prot_alphabet}-k{ksize}.search-contain-matches.sig"),
        txt = os.path.join(out_dir, "search-containment", "{sample}.x.{db_name}.{prot_alphabet}-k{ksize}.search-contain-matches.txt"),
    params:
        alpha_cmd = "", #"--protein", #lambda w: "--" + config["alphabet"], # one alpha at at time. rn using ksize instead
        scaled = config["gather_scaled"],
        ksize = lambda w: int(w.ksize)*3
    resources:
        mem_mb=lambda wildcards, attempt: attempt *20000,
        runtime=6000,
    log: os.path.join(logs_dir, "search-contain", "{sample}.x.{db_name}.{prot_alphabet}-k{ksize}.search-contain.log")
    benchmark: os.path.join(benchmarks_dir, "search-contain", "{sample}.x.{db_name}.{prot_alphabet}-k{ksize}.search-contain.benchmark")
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

rule sourmash_search_containment_nucleotide:
    input:
        nucl_query=rules.sourmash_sketch_nucleotide.output,
        database=os.path.join(database_dir, "{db_name}.{nucl_alphabet}-k{ksize}.sbt.zip")
    output:
        csv = os.path.join(out_dir, "search-containment", "{sample}.x.{db_name}.{nucl_alphabet}-k{ksize}.search-contain-matches.csv"),
        matches = os.path.join(out_dir, "search-containment", "{sample}.x.{db_name}.{nucl_alphabet}-k{ksize}.search-contain-matches.sig"),
        txt = os.path.join(out_dir, "search-containment", "{sample}.x.{db_name}.{nucl_alphabet}-k{ksize}.search-contain-matches.txt"),
    params:
        alpha_cmd = "", #"--protein", #lambda w: "--" + config["alphabet"], # one alpha at at time. rn using ksize instead
        scaled = config["gather_scaled"],
    resources:
        mem_mb=lambda wildcards, attempt: attempt *20000,
        runtime=6000,
    log: os.path.join(logs_dir, "search-contain", "{sample}.x.{db_name}.{nucl_alphabet}-k{ksize}.search-contain.log")
    benchmark: os.path.join(benchmarks_dir, "search-contain", "{sample}.x.{db_name}.{nucl_alphabet}-k{ksize}.search-contain.benchmark")
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

## indexing rules (just working here --> split to separate snakefile)
rule protein_signames_to_file:
    input:  expand(os.path.join(out_dir, "signatures", "{sample}.protein.sig"), sample=sample_names),
    output: os.path.join(out_dir, "index", "{basename}.protein-signatures.txt")
    log: os.path.join(logs_dir, "index", "{basename}.build-siglist-for-index.protein.txt")
    shell:
        """
        for item in {input}
          do
            echo $item >> {output} 2> {log}
          done
        """

rule nucleotide_signames_to_file:
    input:  expand(os.path.join(out_dir, "signatures", "{sample}.nucleotide.sig"), sample=sample_names),
    output: os.path.join(out_dir, "index", "{basename}.nucleotide-signatures.txt")
    log: os.path.join(logs_dir, "index", "{basename}.build-siglist-for-index.nucleotide.txt")
    shell:
        """
        for item in {input}
          do
            echo $item >> {output} 2> {log}
          done
        """

rule sourmash_index_sbt_protein:
    input: os.path.join(out_dir, "index", "{basename}.protein-signatures.txt")
    output: os.path.join(out_dir, "index", "{basename}.{alphabet}-k{ksize}.sbt.zip"),
    threads: 1
    params:
        alpha_cmd = "", # better to use this?
        scaled = config["gather_scaled"],
        ksize = lambda w: int(w.ksize)*3
    resources:
        mem_mb=lambda wildcards, attempt: attempt *50000,
        runtime=6000,
    log: os.path.join(logs_dir, "index", "{basename}.{alphabet}-k{ksize}.index-sbt.log")
    benchmark: os.path.join(benchmarks_dir, "index", "{basename}.{alphabet}-k{ksize}.index-sbt.benchmark")
    conda: "envs/sourmash-dev.yml"
    shell:
        """
        sourmash index {output} --ksize {params.ksize} \
        --scaled {wildcards.scaled} {params.alpha_cmd}  \
        --from-file {input}  2> {log}
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
    resources:
        mem_mb=lambda wildcards, attempt: attempt *10000,
        runtime=6000,
    log: os.path.join(logs_dir, "classify", "{filename}.x.{db_name}.{prot_alphabet}-k{ksize}.guess-tax.log")
    benchmark: os.path.join(benchmarks_dir, "classify", "{filename}.x.{db_name}.{prot_alphabet}-k{ksize}.guess-tax.benchmark")
    conda: "envs/sourmash-dev.yml"
    shell: """
        python -m thumper.guess_taxonomy \
            --genome {input.sample_file} --lineages_csv {input.db_info} \
            --matches_sig {input.matches} \
            --clean {output.clean} --dirty {output.dirty} \
            --report {output.report} --summary {output.csv} \
            --match-rank {params.match_rank} \
            --contig-report {output.contig_report}
    """
