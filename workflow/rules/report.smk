rule vcf_report:
    threads: 8
    input:
        ref="resources/genome.fasta",
        bams=get_batch_bams,
        bcfs=lambda w: expand("results/merged-calls/{group}.{{event}}.fdr-controlled.bcf", group=get_report_batch(w))
    output:
        report(directory("results/vcf-report/{batch}.{event}/"), htmlindex="index.html", caption="../report/calls.rst", category="Variant calls")
    params:
        bcfs=lambda w, input: get_batch_bcfs(w, input),
        bams=lambda w: get_batch_bams(w, True),
        format_field = "DP AF OBS",
        max_read_depth = config["report"]["max_read_depth"],
        js_files="{math} {template}".format(
            math=get_resource("math.min.js"),
            template=get_resource("custom-table-report.js"),
        ),
    log:
        "logs/igv-report/{batch}.{event}.log"
    conda:
        "../envs/rbt.yaml"
    threads: 8
    shell:
        "rbt vcf-report {input.ref} --bams {params.bams} --vcfs {params.bcfs} "
        "--formats {params.format_field} --threads {threads} --infos PROB_*  "
        "-d {params.max_read_depth} --custom-js-files {params.js_files} -- {output}"

