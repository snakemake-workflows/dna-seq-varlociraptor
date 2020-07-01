rule igv_report:
    input:
        bcf="results/merged-calls/{group}.{event}.fdr-controlled.bcf",
        ref="resources/genome.fasta",
        ref_idx="resources/genome.fasta.fai",
        bams=lambda w: get_group_bams(w),
        bais=lambda w: get_group_bams(w, bai=True)
    output:
        report("results/igv-report/{group}.{event}.html", caption="../report/calls.rst", category="Variant calls")
    log:
        "logs/igv-report/{group}.{event}.log"
    params:
        os.path.join(os.path.dirname(workflow.snakefile), "resources/igv-report-template.html")
    conda:
        "../envs/igv-reports.yaml"
    shell:
        "bcftools view {input.bcf} > {output}.vcf;"
        "create_report {output}.vcf {input.ref} --flanking 100 "
        "--info-columns ANN dgiDB_drugs cosmic_LEGACY_ID --info-columns-prefixes PROB_ --sample-columns DP AF OBS"
        " --template {params} --tracks {input.bams} --output {output} --standalone 2>&1 > {log}; "
        "rm {output}.vcf"
