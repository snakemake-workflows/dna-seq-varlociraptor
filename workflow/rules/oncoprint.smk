rule plot_oncoprint:
    input:
        bcfs=get_oncoprint_batch
    output:
        report("results/plots/oncoprint/{batch}.{event}.html", category="Oncoprint", caption="../report/oncoprint.rst")
    params: 
        input=lambda wc, input: expand("{group}={bcf}", group=get_oncoprint_groups, bcf=input.bcfs)
    log:
        "logs/oncoprint/{batch}.{event}.plot.log"
    conda:
        "../envs/rbt.yaml"
    shell:
        "rbt oncoprint {params.input} > {output} &> {log}"
    
