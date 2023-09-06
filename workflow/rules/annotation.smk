rule annotate_candidate_variants:
    input:
        calls="results/candidate-calls/{group}.{caller}.{scatteritem}.bcf",
        cache="resources/vep/cache",
        plugins="resources/vep/plugins",
        fasta=genome,
        fai=genome_fai,
    output:
        calls="results/candidate-calls/{group}.{caller}.{scatteritem}.annotated.bcf",
        stats="results/candidate-calls/{group}.{caller}.{scatteritem}.stats.html",
    params:
        plugins=config["annotations"]["vep"]["candidate_calls"]["plugins"],
        extra="{} --vcf_info_field ANN ".format(
            config["annotations"]["vep"]["candidate_calls"]["params"]
        ),
    log:
        "logs/vep/{group}.{caller}.{scatteritem}.annotate_candidates.log",
    benchmark:
        "benchmarks/vep/{group}.{caller}.{scatteritem}.annotate_candidates.tsv"
    threads: get_vep_threads()
    wrapper:
        "v2.5.0/bio/vep/annotate"


rule annotate_variants:
    input:
        calls="results/calls/{group}.{analysis}.{scatteritem}.bcf",
        cache="resources/vep/cache",
        plugins="resources/vep/plugins",
        revel=lambda wc: get_plugin_aux("REVEL"),
        revel_tbi=lambda wc: get_plugin_aux("REVEL", True),
        fasta=genome,
        fai=genome_fai,
    output:
        calls="results/calls/{group}.{analysis}.{scatteritem}.annotated.bcf",
        stats="results/calls/{group}.{analysis}.{scatteritem}.stats.html",
    params:
        # Pass a list of plugins to use, see https://www.ensembl.org/info/docs/tools/vep/script/vep_plugins.html
        # Plugin args can be added as well, e.g. via an entry "MyPlugin,1,FOO", see docs.
        plugins=config["annotations"]["vep"]["final_calls"]["plugins"],
        extra="{} --vcf_info_field ANN --hgvsg".format(
            config["annotations"]["vep"]["final_calls"]["params"]
        ),
    log:
        "logs/vep/{group}.{analysis}.{scatteritem}.annotate.log",
    threads: get_vep_threads()
    wrapper:
        "v2.5.0/bio/vep/annotate"


# TODO What about multiple ID Fields?
rule annotate_vcfs:
    input:
        bcf="results/calls/{prefix}.bcf",
        annotations=get_annotation_vcfs(),
        idx=get_annotation_vcfs(idx=True),
    output:
        "results/calls/{prefix}.db-annotated.bcf",
    log:
        "logs/annotate-vcfs/{prefix}.log",
    params:
        pipes=get_annotation_pipes,
    conda:
        "../envs/snpsift.yaml"
    threads: 4
    shell:
        "(bcftools view --threads {threads} {input.bcf} {params.pipes} | bcftools view --threads {threads} -Ob > {output}) 2> {log}"


rule annotate_dgidb:
    input:
        "results/calls/{prefix}.bcf",
    params:
        datasources=get_dgidb_datasources(),
    output:
        "results/calls/{prefix}.dgidb.bcf",
    log:
        "logs/annotate-dgidb/{prefix}.log",
    conda:
        "../envs/rbt.yaml"
    resources:
        dgidb_requests=1,
    shell:
        "rbt vcf-annotate-dgidb {input} {params.datasources} > {output} 2> {log}"


rule gather_annotated_calls:
    input:
        calls=get_gather_annotated_calls_input(),
        idx=get_gather_annotated_calls_input(ext="bcf.csi"),
    output:
        "results/final-calls/{group}.{analysis}.annotated.bcf",
    log:
        "logs/gather-annotated-calls/{group}.{analysis}.log",
    params:
        extra="-a",
    wrapper:
        "v2.3.2/bio/bcftools/concat"
