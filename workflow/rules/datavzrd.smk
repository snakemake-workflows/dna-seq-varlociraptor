rule split_call_tables:
    input:
        "results/tables/{group}.{event}.fdr-controlled.tsv",
    output:
        coding="results/tables/{group}.{event}.coding.fdr-controlled.tsv",
        noncoding="results/tables/{group}.{event}.noncoding.fdr-controlled.tsv",
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
        calls=get_oncoprint_input,
        group_annotation=config.get("groups", []),
    output:
        gene_oncoprint="results/tables/oncoprints/{batch}.{event}/gene-oncoprint.tsv",
        gene_oncoprint_sortings=directory(
            "results/tables/oncoprints/{batch}.{event}/label_sortings/"
        ),
        variant_oncoprints=directory(
            "results/tables/oncoprints/{batch}.{event}/variant-oncoprints"
        ),
    log:
        "logs/prepare_oncoprint/{batch}.{event}.log",
    params:
        groups=get_report_batch,
        labels=get_heterogeneous_labels(),
    conda:
        "../envs/oncoprint.yaml"
    script:
        "../scripts/oncoprint.py"


rule render_datavzrd_config:
    input:
        template=workflow.source_path(
            "../resources/datavzrd/variant-calls-template.datavzrd.yaml"
        ),
        variant_oncoprints=get_oncoprint("variant"),
    output:
        "resources/datavzrd/{batch}.{event}.datavzrd.yaml",
    params:
        gene_oncoprint=get_oncoprint("gene"),
        variant_oncoprints=get_variant_oncoprint_tables,
        groups=get_report_batch,
        coding_calls=get_datavzrd_data(impact="coding"),
        noncoding_calls=get_datavzrd_data(impact="noncoding"),
        spec_observations=workflow.source_path(
            "../resources/datavzrd/spec_observations.json"
        ),
        data_observations=workflow.source_path(
            "../resources/datavzrd/data_observations.js"
        ),
        varsome_url=get_varsome_url(),
        samples=samples,
        group_annotations=group_annotation,
        labels=get_heterogeneous_labels(),
        oncoprint_sorted_datasets="results/tables/oncoprints/{batch}.{event}/label_sortings/",
    log:
        "logs/datavzrd_render/{batch}.{event}.log",
    template_engine:
        "yte"


rule datavzrd_variants_calls:
    input:
        coding_calls=get_datavzrd_data(impact="coding"),
        noncoding_calls=get_datavzrd_data(impact="noncoding"),
        spec_observations=workflow.source_path(
            "../resources/datavzrd/spec_observations.json"
        ),
        data_observations=workflow.source_path(
            "../resources/datavzrd/data_observations.js"
        ),
        config="resources/datavzrd/{batch}.{event}.datavzrd.yaml",
        gene_oncoprint=get_oncoprint("gene"),
        variant_oncoprints=get_oncoprint("variant"),
    output:
        report(
            directory("results/datavzrd-report/{batch}.{event}.fdr-controlled"),
            htmlindex="index.html",
            caption="../report/calls.rst",
            category="Variant calls",
            labels=get_datavzrd_report_labels,
            subcategory=get_datavzrd_report_subcategory,
        ),
    log:
        "logs/datavzrd_report/{batch}.{event}.log",
    wrapper:
        "v1.21.1/utils/datavzrd"
