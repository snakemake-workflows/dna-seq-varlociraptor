rule igv_report:
    input:
        bcf="merged-calls/{group}.{event}.fdr-controlled.bcf",
        ref="refs/genome.fasta",
        ref_idx="refs/genome.fasta.fai",
        bams=get_group_bams
    output:
        report("igv-report/{group}.{event}.html", caption="../report/calls.rst", category="Variant calls")
    params:
        os.path.join(os.path.dirname(workflow.snakefile), "resources/igv-report-template.html")
    conda:
        "../envs/igv-reports.yaml"
    shell:
        "bcftools view {input.bcf} > {output}.vcf;"
        "create_report {output}.vcf {input.ref} --flanking 100 "
        "--info-columns ANN dgiDB_drugs cosmic_LEGACY_ID --info-columns-prefixes PROB_ dbNSFP_ --sample-columns DP AF OBS"
        " --template {params} --tracks {input.bams} --output {output} --standalone; "
        "rm {output}.vcf"
