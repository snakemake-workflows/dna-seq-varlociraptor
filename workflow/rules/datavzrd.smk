rule split_call_tables:
    input:
        "results/tables/{group}.{event}.fdr-controlled.tsv",
    output:
        coding="results/tables/{group}.{event}.coding.fdr-controlled.tsv",
        noncoding="results/tables/{group}.{event}.noncoding.fdr-controlled.tsv",
    script:
        "../scripts/split-call-tables.py"


rule render_datavzrd_config:
    input:
        workflow.source_path(
            "../resources/datavzrd/variant-calls-template.datavzrd.yaml"
        ),
    output:
        "resources/datavzrd/all.{event}.datavzrd.yaml",
    params:
        groups=groups,
        coding_calls=get_call_tables("coding"),
        noncoding_calls=get_call_tables("noncoding"),
    template_engine:
        "yte"


rule datavzrd_variants_calls:
    input:
        coding_calls=get_call_tables("coding"),
        noncoding_calls=get_call_tables("noncoding"),
        config="resources/datavzrd/all.{event}.datavzrd.yaml",
    output:
        directory("results/datavzrd-report/all.{event}.fdr-controlled"),
    conda:
        "../envs/datavzrd.yaml"
    shell:
        "datavzrd {input.config} --output {output}"
