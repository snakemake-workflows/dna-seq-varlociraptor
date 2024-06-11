
rule create_mutational_context_file:
    input:
        bcf="results/final-calls/{group}.{event}.variants.fdr-controlled.bcf",
        ref=genome,
        fai=genome_fai,
    output:
        "results/mutational_signatures/context/{group}.{event}.tsv",
    log:
        "logs/mutational_signatures/context/{group}.{event}.log",
    conda:
        "../envs/mutational_context.yaml"
    script:
        "../scripts/create_mutational_context.py"

rule download_cosmic_signatures:
    output:
        "resources/cosmic_signatures.txt",
    params:
        url=lambda wc: "https://cog.sanger.ac.uk/cosmic-signatures-production/documents/COSMIC_v3.4_SBS_{}.txt".format(config["ref"]["build"]),
    log:
        "logs/mutational_signatures/download_cosmic.log",
    conda:
        "../envs/curl.yaml"
    envs:
        "curl {url} -o {output} &> {log}"

rule annotate_mutational_signatures:
    input:
        cosmic_signatures="resources/cosmic_signatures.txt",
        context="results/mutational_signatures/context/{group}.{event}.tsv",
    output:
        "results/mutational_signatures/{group}.{event}.tsv",
    params:
        build=config["ref"]["build"],
    log:
        "logs/mutational_signatures/annotate/{group}.{event}.log",
    conda:
        "../envs/siglasso.yaml"
    script:
        "../scripts/annotate_mutational_signatures.R"