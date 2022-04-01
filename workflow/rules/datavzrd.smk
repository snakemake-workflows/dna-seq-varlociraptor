rule split_call_tables:
    input:
        "results/tables/{group}.{event}.fdr-controlled.tsv",
    output:
        coding="results/tables/{group}.{event}.coding.fdr-controlled.tsv",
        noncoding="results/tables/{group}.{event}.noncoding.fdr-controlled.tsv",
        coding_plot_data="results/tables/{group}.{event}.plotdata.coding.fdr-controlled.tsv",
        noncoding_plot_data="results/tables/{group}.{event}.plotdata.noncoding.fdr-controlled.tsv",
        plot_spec="results/specs/{group}.{event}.varplot.json",
    params:
        sorting=lambda wc: config["calling"]["fdr-control"]["events"][wc.event].get(
            "sort", list()
        ),
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
        "resources/datavzrd/{batch}.{event}.datavzrd.yaml",
    params:
        groups=get_report_batch,
        coding_calls=get_datavzrd_data(impact="coding"),
        noncoding_calls=get_datavzrd_data(impact="noncoding"),
        coding_plotdata=get_datavzrd_data(impact="coding", kind="plotdata"),
        noncoding_plotdata=get_datavzrd_data(impact="noncoding", kind="plotdata"),
        plot_spec=get_datavzrd_data(kind="plotspec"),
        spec_observations=workflow.source_path(
            "../resources/datavzrd/spec_observations.json"
        ),
        data_observations=workflow.source_path(
            "../resources/datavzrd/data_observations.js"
        ),
        varsome_url=get_varsome_url(),
    log:
        "logs/datavzrd_render/{batch}.{event}.log",
    template_engine:
        "yte"


rule datavzrd_variants_calls:
    input:
        coding_calls=get_datavzrd_data(impact="coding"),
        noncoding_calls=get_datavzrd_data(impact="noncoding"),
        coding_plotdata=get_datavzrd_data(impact="coding", kind="plotdata"),
        noncoding_plotdata=get_datavzrd_data(impact="noncoding", kind="plotdata"),
        plot_spec=get_datavzrd_data(kind="plotspec"),
        spec_observations=workflow.source_path(
            "../resources/datavzrd/spec_observations.json"
        ),
        data_observations=workflow.source_path(
            "../resources/datavzrd/data_observations.js"
        ),
        config="resources/datavzrd/{batch}.{event}.datavzrd.yaml",
    output:
        report(
            directory("results/datavzrd-report/{batch}.{event}.fdr-controlled"),
            htmlindex="index.html",
            caption="../report/calls.rst",
            category="Variant calls",
            labels=get_datavzrd_report_labels,
        ),
    conda:
        "../envs/datavzrd.yaml"
    log:
        "logs/datavzrd_report/{batch}.{event}.log",
    shell:
        "datavzrd {input.config} --output {output} &> {log}"
