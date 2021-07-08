## indexing rules 

rule index:
    input:
        expand(os.path.join(out_dir, "index", "{index}.sbt.zip"), index = tp.build_index_names(config))


rule index_sbt:
    input: os.path.join(out_dir, "signatures", "{basename}.signatures.txt")
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
