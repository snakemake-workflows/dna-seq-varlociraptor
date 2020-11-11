rule vcf_report:
    input:
        ref="resources/genome.fasta",
        bams=get_batch_bams,
        bcfs=lambda w: expand("results/merged-calls/{group}.{{event}}.fdr-controlled.bcf", group=get_report_batch(w))
    output:
        report(directory("results/vcf-report/{batch}.{event}/"), htmlindex="index.html", caption="../report/calls.rst", category="Variant calls")
    params:
        bcfs=lambda w: expand("{group}=results/merged-calls/{group}.{event}.fdr-controlled.bcf", group=get_report_batch(w), event=w.event),
        bams=lambda w: get_batch_bams(w, True),
        format_field = "DP AF OBS",
        template = Path(workflow.snakefile).parent / "resources/custom-table-report.js",
        max_read_depth = config["report"]["max_read_depth"],
        js_files = Path(workflow.snakefile).parent / "resources/math.min.js"
    log:
        "logs/igv-report/{batch}.{event}.log"
    conda:
        "../envs/rbt.yaml"
    shell:
        "rbt vcf-report {input.ref} --bams {params.bams} --vcfs {params.bcfs} --format {params.format_field} "
        "--info PROB_* --js {params.template} -d {params.max_read_depth} --js-file {params.js_files} -- {output}"
