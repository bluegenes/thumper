# For later? Expand sample, alpha, ksize, do NOT expand dbname. Aggregate reports over alpha/ksize
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
