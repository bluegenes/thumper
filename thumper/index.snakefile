## indexing rules 

localrules: signames_to_file

rule signames_to_file:
    input:  expand(os.path.join(out_dir, "signatures", "{sample}.sig"), sample=sample_names),
    output: os.path.join(out_dir, "index", "{basename}.signatures.txt")
    run:
        with open(str(output), "w") as outF:
            for inF in input:
                outF.write(str(inF) + "\n")

rule index_sbt:
    input: os.path.join(out_dir, "index", "{basename}.signatures.txt")
    output: os.path.join(out_dir, "index", "{basename}.{alphabet}-k{ksize}-scaled{scaled}.sbt.zip"),
    threads: 1
    params:
        alpha_cmd = lambda w: alphabet_info[w.alphabet]["alpha_cmd"],
        ksize = lambda w: int(w.ksize)*int(alphabet_info[w.alphabet]["ksize_multiplier"]),
    resources:
        mem_mb=lambda wildcards, attempt: attempt *50000,
        runtime=6000,
    log: os.path.join(logs_dir, "index", "{basename}.{alphabet}-k{ksize}-scaled{scaled}.index-sbt.log")
    benchmark: os.path.join(benchmarks_dir, "index", "{basename}.{alphabet}-k{ksize}-scaled{scaled}.index-sbt.benchmark")
    conda: "envs/sourmash-dev.yml"
    shell:
        """
        sourmash index {output} --ksize {params.ksize} \
        --scaled {wildcards.scaled} {params.alpha_cmd}  \
        --from-file {input}  2> {log}
        """

rule index_lca:
    input: 
        sigfile= os.path.join(out_dir, "index", "{basename}.signatures.txt"),
        lineages= config.get('lineages_csv', "")
    output:
        os.path.join(out_dir,"index", "{basename}.{alphabet}-k{ksize}-scaled{scaled}.lca.json.gz"),
    threads: 1
    params:
        alpha_cmd = lambda w: alphabet_info[w.alphabet]["alpha_cmd"],
        ksize = lambda w: int(w.ksize)*int(alphabet_info[w.alphabet]["ksize_multiplier"]),
        report = lambda w: os.path.join(out_dir,"index", f"{w.basename}.{w.alphabet}-k{w.ksize}-scaled{w.scaled}.lca.report"),
    resources:
        mem_mb= lambda wildcards, attempt: attempt *100000,
        runtime=6000,
    log: os.path.join(logs_dir, "index", "{basename}.{alphabet}-k{ksize}-scaled{scaled}.index-lca.log")
    benchmark: os.path.join(benchmarks_dir, "index", "{basename}.{alphabet}-k{ksize}-scaled{scaled}.index-lca.benchmark")
    conda: "envs/sourmash-dev.yml"
    shell: # start column 3 or 4??
        """
        sourmash lca index --ksize {params.ksize} \
          --scaled {wildcards.scaled} --from-file {input.sigfile} \
          --start-column 4 --require-taxonomy --report {params.report} \
          {params.alpha_cmd} {input.lineages} {output} 2> {log} \
        """

