rule split_call_tables:
    input:
        "results/tables/{group}.{event}.fdr-controlled.tsv",
    output:
        coding="results/tables/{group}.{event}.coding.fdr-controlled.tsv",
        noncoding="results/tables/{group}.{event}.noncoding.fdr-controlled.tsv",
        coding_plot="results/tables/{group}.{event}.plotdata.coding.fdr-controlled.tsv",
        noncoding_plot="results/tables/{group}.{event}.plotdata.noncoding.fdr-controlled.tsv",
        plot_spec="results/specs/{group}.{event}.varplot.json",
    log:
        "logs/split_tables/{group}.{event}.log",
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
        coding_calls=get_datavzrd_data(impact="coding"),
        noncoding_calls=get_datavzrd_data(impact="noncoding"),
        coding_plotdata=get_datavzrd_data(impact="coding", kind="plotdata"),
        noncoding_plotdata=get_datavzrd_data(impact="noncoding", kind="plotdata"),
        plot_spec=get_datavzrd_data(kind="plotspec"),
    log:
        "logs/datavzrd_render/{event}.log",
    template_engine:
        "yte"


rule datavzrd_variants_calls:
    input:
        coding_calls=get_datavzrd_data(impact="coding"),
        noncoding_calls=get_datavzrd_data(impact="noncoding"),
        config="resources/datavzrd/all.{event}.datavzrd.yaml",
    output:
        directory("results/datavzrd-report/all.{event}.fdr-controlled"),
    conda:
        "../envs/datavzrd.yaml"
    log:
        "logs/datavzrd_report/{event}.log",
    shell:
        "datavzrd {input.config} --output {output} &> {log}"
