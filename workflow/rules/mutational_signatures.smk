
rule download_mutational_signatures_reference:
    output:
        directory("resources/mutational_signatures_references/tsb")
    params:
        ref_path=lambda wc, output: output[0].rsplit("/", 1)[0],
        build=config["ref"]["build"],
    conda:
        "../envs/sigprofilerassignment.yaml"
    script:
        "../scripts/download_mutational_signatures_reference.py"


rule annotate_mutational_signatures:
    input:
        "results/final-calls/{group}.{event}.variants.fdr-controlled.bcf",
    output:
        directory("results/mutational_signatures/{group}.{event}")
    params:
        build=config["ref"]["build"],
    log:
        "logs/mutational_signatures/annotate/{group}.{event}.log"
    conda:
        "../envs/sigprofilerassignment.yaml"
    script:
        "../scripts/annotate_mutational_signatures.py"