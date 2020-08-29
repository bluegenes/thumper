## indexing rules 
rule signames_to_file:
    input:  expand(os.path.join(out_dir, "signatures", "{sample}.sig"), sample=sample_names),
    output: os.path.join(out_dir, "index", "{basename}.signatures.txt")
    log: os.path.join(logs_dir, "index", "{basename}.build-siglist-for-index.txt")
    shell:
        """
        for item in {input}
          do
            echo $item >> {output} 2> {log}
          done
        """

rule index_sbt:
    input: os.path.join(out_dir, "index", "{basename}.signatures.txt")
    output: os.path.join(out_dir, "index", "{basename}.{alphabet}-k{ksize}.sbt.zip"),
    threads: 1
    params:
        alpha_cmd = lambda w: alphabet_info[w.alphabet]["alpha_cmd"],
        scaled = lambda w: alphabet_info[w.alphabet]["scaled"],
        ksize = lambda w: int(w.ksize)*int(alphabet_info[w.alphabet]["ksize_multiplier"]),
    resources:
        mem_mb=lambda wildcards, attempt: attempt *50000,
        runtime=6000,
    log: os.path.join(logs_dir, "index", "{basename}.{alphabet}-k{ksize}.index-sbt.log")
    benchmark: os.path.join(benchmarks_dir, "index", "{basename}.{alphabet}-k{ksize}.index-sbt.benchmark")
    conda: "envs/sourmash-dev.yml"
    shell:
        """
        sourmash index {output} --ksize {params.ksize} \
        --scaled {params.scaled} {params.alpha_cmd}  \
        --from-file {input}  2> {log}
        """
