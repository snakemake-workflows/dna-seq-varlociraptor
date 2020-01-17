rule vg2svg:
    input:
        "{prefix}.vl.json"
    output:
        "{prefix}.svg"
    conda:
        "../envs/vega.yaml"
    shell:
        "vl2svg {input} {output}"
