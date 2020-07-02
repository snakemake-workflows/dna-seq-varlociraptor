rule plot_oncoprint:
    input:
        bcfs=lambda w: expand("results/merged-calls/{group}.{{event}}.fdr-controlled.bcf", group=get_oncoprint_batch(w))
    output:
        report("results/plots/oncoprint/{batch}.{event}.html", category="Oncoprint", caption="../report/oncoprint.rst")
    params: 
        groups=lambda w, input: expand("{group}={path}", zip, group=get_oncoprint_batch(w), path=input)
    log:
        "logs/oncoprint/{batch}.{event}.plot.log"
    conda:
        "../envs/oncoprint.yaml"
    shell:
        "rbt oncoprint --vep-annotation {params.groups} > {output} 2> {log}"
