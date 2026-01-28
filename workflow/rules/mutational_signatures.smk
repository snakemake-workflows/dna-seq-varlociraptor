rule create_mutational_context_file:
    input:
        bcf="results/final-calls/{group}/{group}.{event}.variants.fdr-controlled.bcf",
        idx="results/final-calls/{group}/{group}.{event}.variants.fdr-controlled.bcf.csi",
        ref=genome,
        fai=genome_fai,
    output:
        context=temp("results/mutational_signatures/{group}.{event}.{sample}.context.tsv"),
        counts=temp("results/mutational_signatures/{group}.{event}.{sample}.counts.tsv"),
    log:
        "logs/mutational_signatures/context/{group}.{event}.{sample}.log",
    params:
        sample_alias=lookup(query="sample_name == '{sample}'", cols="alias", within=samples)
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
        context="results/mutational_signatures/{group}.{event}.{sample}.context.tsv",
    output:
        temp(
            expand(
                "results/mutational_signatures/{{group}}.{{event}}.{{sample}}.{vaf}.tsv",
                vaf=mutational_signature_vaf_thresholds,
            )
        ),
    params:
        build=config["ref"]["build"],
    log:
        "logs/mutational_signatures/annotate/{group}.{event}.{sample}.log",
    conda:
        "../envs/siglasso.yaml"
    script:
        "../scripts/annotate_mutational_signatures.R"


rule join_mutational_signatures:
    input:
        expand(
            "results/mutational_signatures/{{group}}.{{event}}.{{sample}}.{vaf}.tsv",
            vaf=mutational_signature_vaf_thresholds,
        ),
    output:
        temp("results/mutational_signatures/{group}.{event}.{sample}.tsv"),
    log:
        "logs/mutational_signatures/join/{group}.{event}.{sample}.log",
    shell:
        """
        cat <(echo "Signature\tFrequency\tMinimum VAF") {input} >> {output} 2> {log}
        """


rule annotate_descriptions:
    input:
        sig="results/mutational_signatures/{group}.{event}.{sample}.tsv",
        desc=workflow.source_path("../resources/cosmic_signature_desc_v3.4.tsv"),
    output:
        "results/mutational_signatures/{group}.{event}.{sample}.annotated.tsv",
    log:
        "logs/mutational_signatures/annotate/{group}.{event}.{sample}.log",
    conda:
        "../envs/pandas.yaml"
    script:
        "../scripts/annotate_descriptions.py"


rule plot_mutational_signatures:
    input:
        signatures="results/mutational_signatures/{group}.{event}.{sample}.annotated.tsv",
        counts="results/mutational_signatures/{group}.{event}.{sample}.counts.tsv",
    output:
        report(
            "results/plots/mutational_signatures/{group}.{event}.{sample}.html",
            category="Mutational Signatures",
            subcategory="{event}",
            labels={"group": "{group}", "sample": "{sample}"},
        ),
    log:
        "logs/mutational_signatures/{group}.{event}.{sample}.log",
    conda:
        "../envs/pystats.yaml"
    script:
        "../scripts/plot_mutational_signatures.py"
