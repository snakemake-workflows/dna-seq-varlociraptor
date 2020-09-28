rule igv_report:
    input:
        ref="resources/genome.fasta",
        bams=lambda w: get_group_bams(w),
        bcfs= get_merge_calls_input(".bcf")
    output:
        report(directory("results/igv-report/{group}.{event}/"), caption="../report/calls.rst", category="Variant calls")
    params:
        bcfs=lambda w: [bcf.format(group=w.group, event=w.event) for bcf in get_merge_calls_input_report(w, ".bcf")],
        bams=lambda w: [bam.format(group=w.group, event=w.event) for bam in get_group_bams_report(w)]
    log:
        "logs/igv-report/{group}.{event}.log"
    conda:
        "../envs/rbt.yaml"
    shell:
        "rbt vcf-report {input.ref} --bams {params.bams} --vcfs {params.bcfs} -- {output}"
