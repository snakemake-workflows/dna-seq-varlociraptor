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
    conda:
        "../envs/pandas.yaml"
    script:
        "../scripts/split-call-tables.py"


rule prepare_oncoprint:
    input:
        get_oncoprint_input,
    output:
        "results/tables/{batch}.oncoprint.tsv",
    log:
        "logs/prepare_oncoprint/{batch}.log",
    conda:
        "../envs/pandas.yaml"
    script:
        "../scripts/oncoprint.py"


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
        spec_observations=workflow.source_path(
            "../resources/datavzrd/spec_observations.json"
        ),
        data_observations=workflow.source_path(
            "../resources/datavzrd/data_observations.js"
        ),
        varsome_url=get_varsome_url(),
    log:
        "logs/datavzrd_render/{event}.log",
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
        config="resources/datavzrd/all.{event}.datavzrd.yaml",
    output:
        directory("results/datavzrd-report/all.{event}.fdr-controlled"),
    conda:
        "../envs/datavzrd.yaml"
    log:
        "logs/datavzrd_report/{event}.log",
    shell:
        "datavzrd {input.config} --output {output} &> {log}"
