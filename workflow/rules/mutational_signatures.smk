rule create_mutational_context_file:
    input:
        bcf="results/final-calls/{group}.{event}.variants.fdr-controlled.bcf",
        ref=genome,
        fai=genome_fai,
    output:
        context=temp("results/mutational_signatures/{group}.{event}.context.tsv"),
        counts=temp("results/mutational_signatures/{group}.{event}.counts.tsv"),
    log:
        "logs/mutational_signatures/context/{group}.{event}.log",
    conda:
        "../envs/pystats.yaml"
    script:
        "../scripts/create_mutational_context.py"


rule download_cosmic_signatures:
    output:
        "resources/cosmic_signatures.txt",
    params:
        # when updating signature version here, also update workflow/resources/cosmic_signature_desc_v3.4.tsv
        url=lambda wc: "https://cog.sanger.ac.uk/cosmic-signatures-production/documents/COSMIC_v3.4_SBS_{}.txt".format(
            config["ref"]["build"]
        ),
    log:
        "logs/mutational_signatures/download_cosmic.log",
    conda:
        "../envs/curl.yaml"
    shell:
        "curl {params.url} -o {output} &> {log}"


rule annotate_mutational_signatures:
    input:
        cosmic_signatures="resources/cosmic_signatures.txt",
        context="results/mutational_signatures/{group}.{event}.context.tsv",
    output:
        temp(
            expand(
                "results/mutational_signatures/{{group}}.{{event}}.{vaf}.tsv",
                vaf=range(0, 101, 10),
            )
        ),
    params:
        build=config["ref"]["build"],
    log:
        "logs/mutational_signatures/annotate/{group}.{event}.log",
    conda:
        "../envs/siglasso.yaml"
    script:
        "../scripts/annotate_mutational_signatures.R"


rule join_mutational_signatures:
    input:
        expand(
            "results/mutational_signatures/{{group}}.{{event}}.{vaf}.tsv",
            vaf=range(0, 101, 10),
        ),
    output:
        temp("results/mutational_signatures/{group}.{event}.tsv"),
    log:
        "logs/mutational_signatures/join/{group}.{event}.log",
    shell:
        """
        cat <(echo "Signature\tFrequency\tMinimum VAF") {input} >> {output} 2> {log}
        """


rule annotate_descriptions:
    input:
        sig="results/mutational_signatures/{group}.{event}.tsv",
        desc=workflow.source_path("../resources/cosmic_signature_desc_v3.4.tsv"),
    output:
        "results/mutational_signatures/{group}.{event}.annotated.tsv",
    log:
        "logs/mutational_signatures/annotate/{group}.{event}.log",
    conda:
        "../envs/pandas.yaml"
    script:
        "../scripts/annotate_descriptions.py"


rule plot_mutational_signatures:
    input:
        signatures="results/mutational_signatures/{group}.{event}.annotated.tsv",
        counts="results/mutational_signatures/{group}.{event}.counts.tsv",
    output:
        report(
            "results/plots/mutational_signatures/{group}.{event}.html",
            category="Mutational Signatures",
            subcategory="{event}",
            labels={"group": "{group}"},
        ),
    log:
        "logs/mutational_signatures/{group}.{event}.log",
    conda:
        "../envs/pystats.yaml"
    script:
        "../scripts/plot_mutational_signatures.py"
