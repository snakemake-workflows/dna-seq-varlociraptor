rule plot_oncoprint:
    input:
        bcfs=lambda w: expand("results/merged-calls/{group}.{{event}}.fdr-controlled.bcf", group=get_oncoprint_batch(w))
    output:
        report("results/plots/oncoprint/{batch}.{event}.html", category="Oncoprint", caption="../report/oncoprint.rst")
    params: 
        input=lambda w: expand("{group}=results/merged-calls/{group}.{event}.fdr-controlled.bcf", group=get_oncoprint_batch(w), event=w.event)
    log:
        "logs/oncoprint/{batch}.{event}.plot.log"
    conda:
        "../envs/rbt.yaml"
    shell:
        "rbt oncoprint {params.input} > {output} 2> {log}"