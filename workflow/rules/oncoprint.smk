rule plot_oncoprint:
    input:
        bcfs=get_oncoprint_batch
    output:
        report("results/plots/oncoprint/{batch}.{event}.html", category="Oncoprint", caption="../report/oncoprint.rst")
    params: 
        input=lambda wc, input: ["{name}={path}".format(name=bcf_path.split("/")[-1].split(".")[0], path=bcf_path) for bcf_path in input.bcfs]
    log:
        "logs/oncoprint/{batch}.{event}.plot.log"
    conda:
        "../envs/rbt.yaml"
    shell:
        "rbt oncoprint {params.input} > {output} 2> {log}"