
rule preprocess_mutational_variants:
    input:
        "results/final-calls/{group}.{event}.variants.fdr-controlled.bcf",
    output:
        "results/mutational_signatures/{group}.{event}.tsv"
    conda:
        "../envs/pysam.yaml"
    log:
        "logs/mutational_signatures/preprocess_mutational_variants/{group}_{event}.log"
    script:
        "../scripts/preprocess_mutational_variants.py"

rule concat_mutational_calls:
    input:
        expand("results/mutational_signatures/{group}.{{event}}.tsv", group=variants_groups)
    output:
        "results/mutational_signatures/all.{event}.tsv"
    log:
        "logs/mutational_signatures/preprocess_mutational_variants/{event}.log"
    shell:
        """
        echo -e "CHROM\tPOS\tREF\tALT\tPID" > {output}
        cat {input} >> {output} 2> {log}
        """

rule annotate_mutational_signatures:
    input:
        "results/mutational_signatures/all.{event}.tsv"
    output:
        "results/mutational_signatures/{event}.tsv"
    log:
     "logs/mutational_signatures/annotate/{events}.log"
    conda:
        "../envs/yapsa.yaml"
    script:
        "../scripts/annotate_mutational_signatures.R"