rule annotate_variants:
    input:
        calls="results/calls/{group}.bcf",
        cache="resources/vep/cache",
        plugins="resources/vep/plugins"
    output:
        calls="results/calls/{group}.annotated.bcf",
        stats=report("results/calls/{group}.stats.html", caption="../report/stats.rst", category="QC")
    params:
        # Pass a list of plugins to use, see https://www.ensembl.org/info/docs/tools/vep/script/vep_plugins.html
        # Plugin args can be added as well, e.g. via an entry "MyPlugin,1,FOO", see docs.
        plugins=config["annotations"]["vep"]["plugins"],
        extra="{} --vcf_info_field ANN".format(config["annotations"]["vep"]["params"])
    log:
        "logs/vep/{group}.annotate.log"
    wrapper:
        "0.59.2/bio/vep/annotate"


# TODO What about multiple ID Fields?
rule annotate_vcfs:
    threads:
        4
    input:
        bcf="results/calls/{prefix}.bcf",
        annotations=get_annotation_vcfs(),
        idx=get_annotation_vcfs(idx=True)
    output:
        "results/calls/{prefix}.db-annotated.bcf"
    log:
        "logs/annotate-vcfs/{prefix}.log"
    params:
        extra="-Xmx4g",
        pipes=get_annotation_pipes
    conda:
        "../envs/snpsift.yaml"
    shell:
        "(bcftools view --threads {threads} {input.bcf} {params.pipes} | bcftools view --threads {threads} -Ob > {output}) 2> {log}"


rule annotate_dgidb:
    threads:
        4
    input:
        "results/calls/{prefix}.bcf"
    output:
        "results/calls/{prefix}.dgidb.bcf"
    log:
        "logs/annotate-dgidb/{prefix}.log"
    conda:
        "../envs/annotate_dgidb.yaml"
    resources:
        dgidb_requests=1
    shell:
        "rbt vcf-annotate-dgidb {input} > {output} 2> {log}"
