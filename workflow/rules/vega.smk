rule vg2svg:
    input:
        "{prefix}.vl.json",
    output:
        "{prefix}.svg",
    log:
        "logs/vg2svg/{prefix}.log",
    conda:
        "../envs/vega.yaml"
    shell:
        "vl2svg {input} {output} 2> {log}"
