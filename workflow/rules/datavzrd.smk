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
        compact_oncoprint="results/tables/oncoprints/{batch}.{event}/compact-oncoprint.tsv",
        group_sortings="results/tables/oncoprints/{batch}.{event}/group_sortings.json",
    log:
        "logs/prepare_oncoprint/{batch}.{event}.log",
    params:
        groups=get_report_batch,
        labels=get_heterogeneous_labels(),
        compact_oncoprint_header_labels=get_compact_oncoprint_header_labels(),
    conda:
        "../envs/oncoprint.yaml"
    script:
        "../scripts/oncoprint.py"


rule render_compact_oncoprint_spec:
    input:
        template=workflow.source_path("../resources/datavzrd/spec_compact_oncoprint.json.j2"),
        group_sortings="results/tables/oncoprints/{batch}.{event}/group_sortings.json",
    output:
        "resources/datavzrd/{batch}.{event}/spec_compact_oncoprint.sort_by_{sort_label}.json",
    params:
        annotation_labels=get_compact_oncoprint_annotation_labels,
        group_sorting=lambda w, input: json.load(open(input.group_sortings))[w.sort_label],
    template_engine:
        "jinja2"


rule render_datavzrd_config:
    input:
        template=workflow.source_path(
            "../resources/datavzrd/variant-calls-template.datavzrd.yaml"
        ),
        variant_oncoprints=get_oncoprint("variant"),
        compact_oncoprint="results/tables/oncoprints/{batch}.{event}/compact-oncoprint.tsv",
        specs_compact_oncoprint=expand(
            "resources/datavzrd/{{batch}}.{{event}}/spec_compact_oncoprint.sort_by_{sort_label}.json",
            sort_label=get_heterogeneous_labels().index
        ),
    output:
        "resources/datavzrd/{batch}.{event}.datavzrd.yaml",
    params:
        gene_oncoprint=get_oncoprint("gene"),
        compact_oncoprint=get_oncoprint("compact"),
        variant_oncoprints=get_variant_oncoprint_tables,
        groups=get_report_batch,
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
        compact_oncoprint=get_oncoprint("compact"),
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
        "v2.1.1/utils/datavzrd"
