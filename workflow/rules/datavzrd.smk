rule split_call_tables:
    input:
        "results/tables/{group}.{event}.variants.fdr-controlled.tsv",
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


rule process_fusion_call_tables:
    input:
        "results/tables/{group}.{event}.fusions.fdr-controlled.tsv",
    output:
        fusions="results/tables/{group}.{event}.fusions.joined.fdr-controlled.tsv",
    log:
        "logs/join_partner/{group}.{event}.log",
    conda:
        "../envs/pandas.yaml"
    script:
        "../scripts/join_fusion_partner.py"


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
        groups=lambda wc: get_report_batch(wc, "variants"),
        labels=get_heterogeneous_labels(),
    conda:
        "../envs/oncoprint.yaml"
    script:
        "../scripts/oncoprint.py"


rule render_datavzrd_variant_config:
    input:
        template=workflow.source_path(
            "../resources/datavzrd/variant-calls-template.datavzrd.yaml"
        ),
        variant_oncoprints=get_oncoprint("variant"),
    output:
        "resources/datavzrd/{batch}.{event}.variants.datavzrd.yaml",
    params:
        gene_oncoprint=get_oncoprint("gene"),
        variant_oncoprints=get_variant_oncoprint_tables,
        groups=lambda wc: get_report_batch(wc, "variants"),
        coding_calls=get_datavzrd_data(impact="coding"),
        noncoding_calls=get_datavzrd_data(impact="noncoding"),
        spec_observations=workflow.source_path(
            "../resources/datavzrd/spec_observations.json"
        ),
        spec_short_observations=workflow.source_path(
            "../resources/datavzrd/spec_short_observations.json"
        ),
        data_observations=workflow.source_path(
            "../resources/datavzrd/data_observations.js"
        ),
        data_short_observations=workflow.source_path(
            "../resources/datavzrd/data_short_observations.js"
        ),
        varsome_url=get_varsome_url(),
        samples=samples,
        group_annotations=group_annotation,
        labels=get_heterogeneous_labels(),
        oncoprint_sorted_datasets="results/tables/oncoprints/{batch}.{event}/label_sortings/",
    log:
        "logs/datavzrd_render/{batch}.{event}.variants.log",
    template_engine:
        "yte"


rule render_datavzrd_fusions_config:
    input:
        template=workflow.source_path(
            "../resources/datavzrd/fusion-calls-template.datavzrd.yaml"
        ),
    output:
        "resources/datavzrd/{batch}.{event}.fusions.datavzrd.yaml",
    params:
        groups=lambda wc: get_report_batch(wc, "fusions"),
        fusion_calls=get_datavzrd_data(impact="fusions"),
        spec_observations=workflow.source_path(
            "../resources/datavzrd/spec_observations.json"
        ),
        spec_short_observations=workflow.source_path(
            "../resources/datavzrd/spec_short_observations.json"
        ),
        data_observations=workflow.source_path(
            "../resources/datavzrd/data_observations.js"
        ),
        data_short_observations=workflow.source_path(
            "../resources/datavzrd/data_short_observations.js"
        ),
        samples=samples,
    log:
        "logs/datavzrd_render/{batch}.{event}.fusions.log",
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
        config="resources/datavzrd/{batch}.{event}.variants.datavzrd.yaml",
        gene_oncoprint=get_oncoprint("gene"),
        variant_oncoprints=get_oncoprint("variant"),
    output:
        report(
            directory(
                "results/datavzrd-report/{batch}.{event}.variants.fdr-controlled"
            ),
            htmlindex="index.html",
            caption="../report/calls.rst",
            category="Variant calls",
            labels=lambda wc: get_datavzrd_report_labels(wc, "variants"),
            subcategory=get_datavzrd_report_subcategory,
        ),
    log:
        "logs/datavzrd_report/{batch}.{event}.log",
    wrapper:
        "v2.6.0/utils/datavzrd"


rule datavzrd_fusion_calls:
    input:
        fusion_calls=get_datavzrd_data(impact="fusions"),
        spec_observations=workflow.source_path(
            "../resources/datavzrd/spec_observations.json"
        ),
        data_observations=workflow.source_path(
            "../resources/datavzrd/data_observations.js"
        ),
        config="resources/datavzrd/{batch}.{event}.fusions.datavzrd.yaml",
    output:
        report(
            directory(
                "results/datavzrd-report/{batch}.{event}.fusions.fdr-controlled"
            ),
            htmlindex="index.html",
            caption="../report/calls.rst",
            category="Fusion calls",
            labels=lambda wc: get_datavzrd_report_labels(wc, "fusions"),
            subcategory=get_datavzrd_report_subcategory,
        ),
    log:
        "logs/datavzrd_report/{batch}.{event}.log",
    wrapper:
        "v2.6.0/utils/datavzrd"
