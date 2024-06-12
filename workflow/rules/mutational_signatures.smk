rule create_mutational_context_file:
    input:
        bcf="results/final-calls/{group}.{event}.variants.fdr-controlled.bcf",
        ref=genome,
        fai=genome_fai,
    output:
        "results/mutational_signatures/context/{group}.{event}.{vaf}.tsv",
    log:
        "logs/mutational_signatures/context/{group}.{event}.{vaf}.log",
    conda:
        "../envs/mutational_context.yaml"
    script:
        "../scripts/create_mutational_context.py"


rule download_cosmic_signatures:
    output:
        "resources/cosmic_signatures.txt",
    params:
        url=lambda wc: "https://cog.sanger.ac.uk/cosmic-signatures-production/documents/COSMIC_v3.4_SBS_{}.txt".format(
            config["ref"]["build"]
        ),
    log:
        "logs/mutational_signatures/download_cosmic.log",
    conda:
        "../envs/curl.yaml"
    shell:
        "curl {url} -o {output} &> {log}"


rule annotate_mutational_signatures:
    input:
        cosmic_signatures="resources/cosmic_signatures.txt",
        context="results/mutational_signatures/context/{group}.{event}.{vaf}.tsv",
    output:
        "results/mutational_signatures/{group}.{event}.{vaf}.tsv",
    params:
        build=config["ref"]["build"],
    log:
        "logs/mutational_signatures/annotate/{group}.{event}.{vaf}.log",
    conda:
        "../envs/siglasso.yaml"
    script:
        "../scripts/annotate_mutational_signatures.R"


rule join_mutational_signatures:
    input:
        expand(
            "results/mutational_signatures/{{group}}.{{event}}.{vaf}.tsv",
            vaf=range(10, 101, 10),
        ),
    output:
        "results/mutational_signatures/{group}.{event}.tsv",
    log:
        "logs/mutational_signatures/join/{group}.{event}.log",
    script:
        "../scripts/join_mutational_signatures.py"


rule plot_mutational_signatures:
    input:
        "results/mutational_signatures/{group}.{event}.tsv",
    output:
        report(
            "results/plots/mutational-signatures/{group}.{event}.mutational-burden.svg",
            category="Mutational Signatures",
            subcategory="{group}",
            labels={"event": "{event}"},
        ),
    log:
        "logs/mutational_signatures/{group}.{event}.log",
    conda:
        "../envs/plot_ms.yaml"
    script:
        "scripts/plot_mutational_signatures.py"
