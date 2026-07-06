rule annotate_candidate_variants:
    input:
        calls="results/candidate-calls/{caller}/{group}/{group}.{scatteritem}.bcf",
        cache=access.random("resources/vep/cache"),
        plugins=access.random("resources/vep/plugins"),
        fasta=access.random(genome),
        fai=genome_fai,
    output:
        calls="results/candidate-calls/{caller}/{group}/{group}.{scatteritem}.annotated.bcf",
        stats="results/candidate-calls/{caller}/{group}/{group}.{scatteritem}.stats.html",
    log:
        "logs/vep/{caller}/{group}/{group}.{scatteritem}.annotate_candidates.log",
    benchmark:
        "benchmarks/vep/{caller}/{group}/{group}.{scatteritem}.annotate_candidates.tsv"
    group:
        "candidate-annotation"
    threads: 4
    params:
        plugins=config["annotations"]["vep"]["candidate_calls"]["plugins"],
        extra="{} --vcf_info_field ANN ".format(
            config["annotations"]["vep"]["candidate_calls"]["params"]
        ),
    wrapper:
        "v8.0.0/bio/vep/annotate"


rule annotate_variants:
    input:
        calls="results/calls/varlociraptor/{group}/{group}.{calling_type}.{scatteritem}.bcf",
        cache=access.random("resources/vep/cache"),
        plugins=access.random("resources/vep/plugins"),
        revel=lambda wc: get_plugin_aux("REVEL"),
        revel_tbi=lambda wc: get_plugin_aux("REVEL", index=True),
        cadd_snv=lambda wc: get_plugin_aux("CADD", cadd_variant_type="snv"),
        cadd_snv_tbi=lambda wc: get_plugin_aux(
            "CADD", cadd_variant_type="snv", index=True
        ),
        cadd_indel=lambda wc: get_plugin_aux("CADD", cadd_variant_type="indels"),
        cadd_indel_tbi=lambda wc: get_plugin_aux(
            "CADD", cadd_variant_type="indels", index=True
        ),
        fasta=access.random(genome),
        fai=genome_fai,
    output:
        calls="results/calls/vep_annotated/{group}/{group}.{calling_type}.{scatteritem}.bcf",
        stats="results/calls/vep_annotated/{group}/{group}.{calling_type}.{scatteritem}.stats.html",
    log:
        "logs/vep/{group}.{calling_type}.{scatteritem}.annotate.log",
    group:
        "annotation"
    threads: 4
    params:
        # Pass a list of plugins to use, see https://www.ensembl.org/info/docs/tools/vep/script/vep_plugins.html
        # Plugin args can be added as well, e.g. via an entry "MyPlugin,1,FOO", see docs.
        plugins=lambda wc: [
            p.replace(
                "CADD",
                f"CADD,snv={get_plugin_aux('CADD', cadd_variant_type= 'snv')},indels={get_plugin_aux('CADD', cadd_variant_type= 'indels')}",
            )
            for p in config["annotations"]["vep"]["final_calls"]["plugins"]
        ],
        extra="{} --vcf_info_field ANN --hgvsg".format(
            config["annotations"]["vep"]["final_calls"]["params"]
        ),
    wrapper:
        "v8.0.0/bio/vep/annotate"


# TODO What about multiple ID Fields?
rule annotate_vcfs:
    input:
        bcf="results/calls/vep_annotated/{prefix}.bcf",
        annotations=get_annotation_vcfs(),
        idx=get_annotation_vcfs(idx=True),
    output:
        "results/calls/db_annotated/{prefix}.bcf",
    log:
        "logs/annotate-vcfs/{prefix}.log",
    group:
        "annotation"
    conda:
        "../envs/snpsift.yaml"
    threads: 4
    params:
        pipes=get_annotation_pipes,
    shell:
        "(bcftools view --threads {threads} {input.bcf} {params.pipes} | "
        "bcftools view --threads {threads} -Ob > {output}) 2> {log}"


rule annotate_dgidb:
    input:
        get_annotate_dgidb_input,
    output:
        "results/calls/dgidb_annotated/{prefix}.bcf",
    log:
        "logs/annotate-dgidb/{prefix}.log",
    group:
        "annotation"
    conda:
        "../envs/rbt.yaml"
    resources:
        dgidb_requests=1,
    params:
        datasources=get_dgidb_datasources(),
    shell:
        "rbt vcf-annotate-dgidb {input} {params.datasources} > {output} 2> {log}"


use rule bcf_index as annotated_index with:
    group:
        "annotation"


ruleorder: annotated_index > bcf_index


rule gather_annotated_calls:
    input:
        calls=get_gather_annotated_calls_input(),
        idx=get_gather_annotated_calls_input(ext="bcf.csi"),
    output:
        "results/final-calls/{group}/{group}.{calling_type}.annotated.bcf",
    log:
        "logs/gather-annotated-calls/{group}/{group}.{calling_type}.log",
    group:
        "annotation"
    params:
        extra="-a",
    wrapper:
        "v2.3.2/bio/bcftools/concat"
