rule igv_report:
    input:
        bcf="results/merged-calls/{group}.{event}.fdr-controlled.bcf",
        ref="results/refs/genome.fasta",
        ref_idx="results/refs/genome.fasta.fai",
        bams=get_group_bams,
        bais=get_group_bais
    output:
        report("results/igv-report/{group}.{event}.html", caption="../results/report/calls.rst", category="Variant calls")
    params:
        os.path.join(os.path.dirname(workflow.snakefile), "resources/igv-report-template.html")
    conda:
        "../envs/igv-reports.yaml"
    log:
        "logs/igv-report/{group}.{event}.log"
    shell:
        "bcftools view {input.bcf} > {output}.vcf;"
        "create_report {output}.vcf {input.ref} --flanking 100 "
        "--info-columns ANN dgiDB_drugs cosmic_LEGACY_ID --info-columns-prefixes PROB_ dbNSFP_ --sample-columns DP AF OBS"
        " --template {params} --tracks {input.bams} --output {output} --standalone 2>&1 > {log}; "
        "rm {output}.vcf"
