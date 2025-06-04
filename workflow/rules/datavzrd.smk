rule split_call_tables:
    input:
        calls="results/tables/{group}.{any_event}.variants.fdr-controlled.tsv",
        population_db=get_cleaned_population_db(),
        population_db_idx=get_cleaned_population_db(idx=True),
    output:
        coding="results/tables/{group}.{any_event}.coding.fdr-controlled.tsv",
        noncoding="results/tables/{group}.{any_event}.noncoding.fdr-controlled.tsv",
    params:
        sorting=lambda wc: lookup(within=get_any_event(wc), dpath="sort", default=list()),
    log:
        "logs/split_tables/{group}.{any_event}.log",
    conda:
        "../envs/split_call_tables.yaml"
    script:
        "../scripts/split-call-tables.py"


rule process_fusion_call_tables:
    input:
        "results/tables/{group}.{any_event}.fusions.fdr-controlled.tsv",
    output:
        fusions="results/tables/{group}.{any_event}.fusions.joined.fdr-controlled.tsv",
    log:
        "logs/join_partner/{group}.{any_event}.log",
    conda:
        "../envs/pandas.yaml"
    script:
        "../scripts/join_fusion_partner.py"


rule prepare_oncoprint:
    input:
        calls=get_oncoprint_input,
        group_annotation=config.get("groups", []),
    output:
        gene_oncoprint="results/tables/oncoprints/{batch}.{any_event}/gene-oncoprint.tsv",
        gene_oncoprint_sortings=directory(
            "results/tables/oncoprints/{batch}.{any_event}/label_sortings/"
        ),
        variant_oncoprints=directory(
            "results/tables/oncoprints/{batch}.{any_event}/variant-oncoprints"
        ),
    log:
        "logs/prepare_oncoprint/{batch}.{any_event}.log",
    params:
        groups=get_report_batch("variants"),
        labels=get_heterogeneous_labels(),
    conda:
        "../envs/oncoprint.yaml"
    script:
        "../scripts/oncoprint.py"


rule datavzrd_variants_calls:
    input:
        coding_calls=get_datavzrd_data(impact="coding"),
        noncoding_calls=get_datavzrd_data(impact="noncoding"),
        linkouts=workflow.source_path("../resources/datavzrd/linkouts.js"),
        spec_observations=workflow.source_path(
            "../resources/datavzrd/spec_observations.json"
        ),
        data_observations=workflow.source_path(
            "../resources/datavzrd/data_observations.js"
        ),
        spec_short_observations=workflow.source_path(
            "../resources/datavzrd/spec_short_observations.json"
        ),
        data_short_observations=workflow.source_path(
            "../resources/datavzrd/data_short_observations.js"
        ),
        config=workflow.source_path(
            "../resources/datavzrd/variant-calls-template.datavzrd.yaml"
        ),
        gene_oncoprint=get_oncoprint("gene"),
        variant_oncoprints=get_oncoprint("variant"),
        oncoprint_sorted_datasets="results/tables/oncoprints/{batch}.{any_event}/label_sortings/",
    output:
        report(
            directory(
                "results/datavzrd-report/{batch}.{any_event}.variants.fdr-controlled"
            ),
            htmlindex="index.html",
            caption="../report/calls.rst",
            category="Variant calls",
            labels=get_datavzrd_report_labels,
            subcategory=get_datavzrd_report_subcategory,
        ),
    log:
        "logs/datavzrd_report/{batch}.{any_event}.log",
    params:
        variant_oncoprints=get_variant_oncoprint_tables,
        groups=get_report_batch("variants"),
        coding_calls=get_datavzrd_data(impact="coding"),
        noncoding_calls=get_datavzrd_data(impact="noncoding"),
        build=config["ref"]["build"],
        samples=samples,
        group_annotations=group_annotation,
        labels=get_heterogeneous_labels(),
        event_desc=lambda wc: lookup(within=get_any_event(wc), dpath="desc", default="No description provided"),
    wrapper:
        "v5.6.1/utils/datavzrd"


rule datavzrd_fusion_calls:
    input:
        fusion_calls=get_datavzrd_data(impact="fusions"),
        spec_observations=workflow.source_path(
            "../resources/datavzrd/spec_observations.json"
        ),
        data_observations=workflow.source_path(
            "../resources/datavzrd/data_observations.js"
        ),
        spec_short_observations=workflow.source_path(
            "../resources/datavzrd/spec_short_observations.json"
        ),
        data_short_observations=workflow.source_path(
            "../resources/datavzrd/data_short_observations.js"
        ),
        config=workflow.source_path(
            "../resources/datavzrd/fusion-calls-template.datavzrd.yaml"
        ),
    output:
        report(
            directory(
                "results/datavzrd-report/{batch}.{any_event}.fusions.fdr-controlled"
            ),
            htmlindex="index.html",
            caption="../report/calls.rst",
            category="Fusion calls",
            labels=get_datavzrd_report_labels,
            subcategory=get_datavzrd_report_subcategory,
        ),
    log:
        "logs/datavzrd_report/{batch}.{any_event}.log",
    params:
        groups=get_report_batch("fusions"),
        samples=samples,
    wrapper:
        "v5.6.1/utils/datavzrd"


rule bedtools_merge:
    input:
        left="results/regions/{group}/{sample}.regions.bed.gz",
        right="results/regions/{group}.covered_regions.bed",
    output:
        "results/coverage/{group}/{sample}.regions.filtered.bed",
    params:
        ## Add optional parameters
        extra="-wa",
    log:
        "logs/bedtools/{group}/{sample}.log",
    wrapper:
        "v2.6.0/bio/bedtools/intersect"


rule coverage_table:
    input:
        lambda wc: expand(
            "results/coverage/{{group}}/{sample}.regions.filtered.bed",
            sample=get_group_samples(wc.group),
        ),
    output:
        "results/coverage/{group}.csv",
    params:
        min_cov=config["gene_coverage"].get("min_avg_coverage", 0),
    conda:
        "../envs/pandas.yaml"
    log:
        "logs/coverage/{group}_coverage_table.log",
    script:
        "../scripts/coverage_table.py"


rule datavzrd_coverage:
    input:
        csv="results/coverage/{group}.csv",
        config=workflow.source_path(
            "../resources/datavzrd/gene-coverage-template.datavzrd.yaml"
        ),
    output:
        report(
            directory("results/datavzrd-report/{group}.coverage"),
            htmlindex="index.html",
            category="Mean read depth per gene",
            labels=lambda wc: {"Group": wc.group},
        ),
    log:
        "logs/datavzrd_report/{group}.coverage.log",
    params:
        samples=lambda wc: get_group_samples(wc.group),
    wrapper:
        "v5.6.1/utils/datavzrd"
