rule vcf_report:
    input:
        ref="resources/genome.fasta",
        bams=lambda w: get_batch_bams(w),
        bcfs=lambda w: expand("results/merged-calls/{group}.{{event}}.fdr-controlled.bcf", group=get_report_batch(w))
    output:
        report(directory("results/vcf-report/{batch}.{event}/"), htmlindex="index.html", caption="../report/calls.rst", category="Variant calls")
    params:
        bcfs=lambda w: expand("{group}=results/merged-calls/{group}.{event}.fdr-controlled.bcf", group=get_report_batch(w), event=w.event),
        bams=lambda w: get_batch_bams(w, True)
    log:
        "logs/igv-report/{batch}.{event}.log"
    conda:
        "../envs/rbt.yaml"
    shell:
        "rbt vcf-report {input.ref} --bams {params.bams} --vcfs {params.bcfs} -- {output}"