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
