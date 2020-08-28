## indexing rules 
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
    output: os.path.join(out_dir, "index", "{basename}.{prot_alphabet}-k{ksize}.sbt.zip"),
    threads: 1
    params:
        alpha_cmd = "", # better to use this?
        scaled = config["protein_scaled"],
        ksize = lambda w: int(w.ksize)*3
    resources:
        mem_mb=lambda wildcards, attempt: attempt *50000,
        runtime=6000,
    log: os.path.join(logs_dir, "index", "{basename}.{prot_alphabet}-k{ksize}.index-sbt.log")
    benchmark: os.path.join(benchmarks_dir, "index", "{basename}.{prot_alphabet}-k{ksize}.index-sbt.benchmark")
    conda: "envs/sourmash-dev.yml"
    shell:
        """
        sourmash index {output} --ksize {params.ksize} \
        --scaled {wildcards.scaled} {params.alpha_cmd}  \
        --from-file {input}  2> {log}
        """

rule sourmash_index_sbt_nucleotide:
    input: os.path.join(out_dir, "index", "{basename}.nucleotide-signatures.txt")
    output: os.path.join(out_dir, "index", "{basename}.{nucl_alphabet}-k{ksize}.sbt.zip"),
    threads: 1
    params:
        alpha_cmd = "", # better to use this?
        scaled = config["nucleotide_scaled"],
        ksize = lambda w: int(w.ksize)
    resources:
        mem_mb=lambda wildcards, attempt: attempt *50000,
        runtime=6000,
    log: os.path.join(logs_dir, "index", "{basename}.{nucl_alphabet}-k{ksize}.index-sbt.log")
    benchmark: os.path.join(benchmarks_dir, "index", "{basename}.{nucl_alphabet}-k{ksize}.index-sbt.benchmark")
    conda: "envs/sourmash-dev.yml"
    shell:
        """
        sourmash index {output} --ksize {params.ksize} \
        --scaled {wildcards.scaled} {params.alpha_cmd}  \
        --from-file {input}  2> {log}
        """
