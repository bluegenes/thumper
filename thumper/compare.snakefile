## compare a group of signatures

localrules: compare_signames_to_file

rule compare_signames_to_file:
    input:  expand(os.path.join(out_dir, "signatures", "{sample}.sig"), sample=sample_names),
    output: os.path.join(out_dir, "compare", "{basename}.signatures.txt")
    run:
        with open(str(output), "w") as outF:
            for inF in input:
                outF.write(str(inF) + "\n")

rule jaccard_compare_sigs:
    input: os.path.join(out_dir, "compare", "{basename}.signatures.txt")
    output: 
        csv=os.path.join(out_dir, "compare", "{basename}.{alphabet}-k{ksize}-scaled{scaled}.jaccard.csv"),
        np=os.path.join(out_dir, "compare", "{basename}.{alphabet}-k{ksize}-scaled{scaled}.jaccard.np"),
    threads: 1
    params:
        alpha_cmd = lambda w: alphabet_info[w.alphabet]["alpha_cmd"],
        ksize = lambda w: int(w.ksize)*int(alphabet_info[w.alphabet]["ksize_multiplier"]),
    resources:
        mem_mb=lambda wildcards, attempt: attempt *10000,
        runtime=6000,
    log: os.path.join(logs_dir, "compare", "{basename}.{alphabet}-k{ksize}-scaled{scaled}.jaccard.log")
    benchmark: os.path.join(benchmarks_dir, "compare", "{basename}.{alphabet}-k{ksize}-scaled{scaled}.jaccard.benchmark")
    conda: "envs/sourmash-dev.yml"
    shell:
        """
        sourmash compare -o {output.np} --csv {output.csv} \
        --ignore-abundance --ksize {params.ksize} \
        {params.alpha_cmd}  \
        --from-file {input}  2> {log}
        """

rule cosine_compare_sigs:
    input: os.path.join(out_dir, "compare", "{basename}.signatures.txt")
    output: 
        csv=os.path.join(out_dir, "compare", "{basename}.{alphabet}-k{ksize}-scaled{scaled}.cosine.csv"),
        np=os.path.join(out_dir, "compare", "{basename}.{alphabet}-k{ksize}-scaled{scaled}.cosine.np"),
    threads: 1
    params:
        alpha_cmd = lambda w: alphabet_info[w.alphabet]["alpha_cmd"],
        ksize = lambda w: int(w.ksize)*int(alphabet_info[w.alphabet]["ksize_multiplier"]),
    resources:
        mem_mb=lambda wildcards, attempt: attempt *10000,
        runtime=6000,
    log: os.path.join(logs_dir, "compare", "{basename}.{alphabet}-k{ksize}-scaled{scaled}.cosine.log")
    benchmark: os.path.join(benchmarks_dir, "compare", "{basename}.{alphabet}-k{ksize}-scaled{scaled}.cosine.benchmark")
    conda: "envs/sourmash-dev.yml"
    shell:
        """
        sourmash compare -o {output.np} --csv {output.csv} \
        --ksize {params.ksize} \
        {params.alpha_cmd}  \
        --from-file {input}  2> {log}
        """

rule containment_compare_sigs:
    input: os.path.join(out_dir, "compare", "{basename}.signatures.txt")
    output: 
        csv=os.path.join(out_dir, "compare", "{basename}.{alphabet}-k{ksize}-scaled{scaled}.containment.csv"),
        np=os.path.join(out_dir, "compare", "{basename}.{alphabet}-k{ksize}-scaled{scaled}.containment.np"),
    threads: 1
    params:
        alpha_cmd = lambda w: alphabet_info[w.alphabet]["alpha_cmd"],
        ksize = lambda w: int(w.ksize)*int(alphabet_info[w.alphabet]["ksize_multiplier"]),
    resources:
        mem_mb=lambda wildcards, attempt: attempt *10000,
        runtime=6000,
    log: os.path.join(logs_dir, "compare", "{basename}.{alphabet}-k{ksize}-scaled{scaled}.containment.log")
    benchmark: os.path.join(benchmarks_dir, "compare", "{basename}.{alphabet}-k{ksize}-scaled{scaled}.containment.benchmark")
    conda: "envs/sourmash-dev.yml"
    shell:
        # containment analyses ignore abundances by default, right? only set coverage..
        """
        sourmash compare -o {output.np} --csv {output.csv} \
        --containment --ignore-abundance --ksize {params.ksize} \
        {params.alpha_cmd}  \
        --from-file {input}  2> {log}
        """

rule sourmash_plot:
    input: os.path.join(out_dir, "compare", "{basename}.{alphabet}-k{ksize}-scaled{scaled}.{comparison}.np"),
    output: os.path.join(out_dir, "compare", "{basename}.{alphabet}-k{ksize}-scaled{scaled}.{comparison}.np.matrix.pdf")
    params:
        plot_dir=lambda w: os.path.join(out_dir, "compare")
    resources:
        mem_mb=lambda wildcards, attempt: attempt *10000,
        runtime=6000,
    conda: "envs/sourmash-dev.yml"
    shell:
        """
        sourmash plot --output-dir {params.plot_dir} --labels --pdf {input}
        """
