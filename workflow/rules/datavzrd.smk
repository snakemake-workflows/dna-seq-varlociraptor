rule process_call_tables:
    input:
        calls="results/tables/{group}/{group}.{event}.variants.fdr-controlled.tsv",
        population_db=get_cleaned_population_db(),
        population_db_idx=get_cleaned_population_db(idx=True),
    output:
        "results/tables/{group}/{group}.{event}.variants.postprocessed.fdr-controlled.tsv",
    log:
        "logs/process_call_tables/{group}.{event}.log",
    conda:
        "../envs/process_call_tables.yaml"
    params:
        sorting=lambda wc: config["calling"]["fdr-control"]["events"][wc.event].get(
            "sort", list()
        ),
    script:
        "../scripts/process-call-tables.py"


rule process_fusion_call_tables:
    input:
        varlociraptor="results/tables/{group}/{group}.{event}.fusions.fdr-controlled.tsv",
        arriba=expand(
            "results/arriba/{sample}.fusions.annotated.tsv",
            sample=lookup(
                within=samples,
                query="group == '{group}' & calling == 'fusions' & datatype == 'rna'",
                cols="sample_name",
            ),
        ),
    output:
        fusions="results/tables/{group}/{group}.{event}.fusions.joined.fdr-controlled.tsv",
    log:
        "logs/join_partner/{group}.{event}.log",
    conda:
        "../envs/pandas.yaml"
    script:
        "../scripts/create_fusions_table_per_group.py"


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
    conda:
        "../envs/oncoprint.yaml"
    params:
        groups=get_report_batch("variants"),
        labels=get_heterogeneous_labels(),
    script:
        "../scripts/oncoprint.py"


rule datavzrd_variants_calls:
    input:
        calls=get_datavzrd_data(calling_type="variants"),
        linkouts=workflow.source_path("../resources/datavzrd/linkouts.js"),
        config=workflow.source_path(
            "../resources/datavzrd/variant-calls-template.datavzrd.yaml"
        ),
        gene_oncoprint=get_oncoprint("gene"),
        variant_oncoprints=get_oncoprint("variant"),
        oncoprint_sorted_datasets="results/tables/oncoprints/{batch}.{event}/label_sortings/",
    output:
        report(
            directory(
                "results/datavzrd-report/{batch}.{event}.variants.fdr-controlled"
            ),
            htmlindex="index.html",
            caption="../report/calls.rst",
            category="Variant calls",
            labels=get_datavzrd_report_labels,
            subcategory=get_datavzrd_report_subcategory,
        ),
    log:
        "logs/datavzrd_report/{batch}.{event}.log",
    params:
        variant_oncoprints=get_variant_oncoprint_tables,
        groups=get_report_batch("variants"),
        build=config["ref"]["build"],
        genebe_genome_build=genebe_genome_build,
        samples=samples,
        group_annotations=group_annotation,
        labels=get_heterogeneous_labels(),
        event_desc=lookup(
            dpath="calling/fdr-control/events/{event}/desc", within=config
        ),
    wrapper:
        "v9.10.1/utils/datavzrd"


rule datavzrd_fusion_calls:
    input:
        fusion_calls=get_datavzrd_data(calling_type="fusions"),
        config=workflow.source_path(
            "../resources/datavzrd/fusion-calls-template.datavzrd.yaml"
        ),
    output:
        report(
            directory(
                "results/datavzrd-report/{batch}.{event}.fusions.fdr-controlled"
            ),
            htmlindex="index.html",
            caption="../report/calls.rst",
            category="Fusion calls",
            labels=get_datavzrd_report_labels,
            subcategory=get_datavzrd_report_subcategory,
        ),
    log:
        "logs/datavzrd_report/{batch}.{event}.log",
    params:
        groups=get_report_batch("fusions"),
        species=lookup(within=config, dpath="ref/species"),
        samples=samples,
    wrapper:
        "v9.10.1/utils/datavzrd"


rule bedtools_merge:
    input:
        left="results/regions/{group}/{sample}.regions.bed.gz",
        right="results/regions/{group}.covered_regions.bed",
    output:
        "results/coverage/{group}/{sample}.regions.filtered.bed",
    log:
        "logs/bedtools/{group}/{sample}.log",
    params:
        ## Add optional parameters
        extra="-wa",
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
    log:
        "logs/coverage/{group}_coverage_table.log",
    conda:
        "../envs/pandas.yaml"
    params:
        min_cov=config["gene_coverage"].get("min_avg_coverage", 0),
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
        "v9.10.1/utils/datavzrd"
