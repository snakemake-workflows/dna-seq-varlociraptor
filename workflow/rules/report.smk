rule vcf_report:
    threads: 8
    input:
        ref="resources/genome.fasta",
        bams=get_batch_bams,
        bcfs=get_batch_bcfs,
    output:
        report(
            directory("results/vcf-report/{batch}.{event}/"),
            htmlindex="index.html",
            caption="../report/calls.rst",
            category="Variant calls",
        ),
    params:
        bcfs=get_report_bcf_params,
        bams=get_report_bam_params,
        format_field="DP AF SOBS OBS AFD",
        max_read_depth=config["report"]["max_read_depth"],
        js_files="{math} {template}".format(
            math=get_resource("math.min.js"),
            template=get_resource("custom-table-report.js"),
        ),
    log:
        "logs/igv-report/{batch}.{event}.log",
    conda:
        "../envs/rbt.yaml"
    threads: 8
    shell:
        "rbt vcf-report {input.ref} --bams {params.bams} --vcfs {params.bcfs} "
        "--formats {params.format_field} --threads {threads} --infos PROB_*  "
        "-d {params.max_read_depth} --custom-js-files {params.js_files} -- {output} 2> {log}"
