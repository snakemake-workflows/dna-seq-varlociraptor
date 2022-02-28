rule vembrane_table:
    input:
        bcf="results/final-calls/{group}.{event}.fdr-controlled.normal-probs.bcf",
    output:
        bcf="results/tables/{group}.{event}.fdr-controlled.tsv",
    conda:
        "../envs/vembrane.yaml"
    params:
        config=get_vembrane_config,
    log:
        "logs/vembrane-table/{group}.{event}.log",
    shell:
        'vembrane table --header "{params.config[header]}" "{params.config[expr]}" {input.bcf} > {output.bcf} 2> {log}'


rule tsv_to_excel:
    input:
        tsv="results/{x}.tsv",
    output:
        xlsx="results/{x}.xlsx",
    conda:
        "../envs/excel.yaml"
    log:
        "logs/tsv_to_xlsx/{x}.log",
    script:
        "../scripts/tsv_to_xlsx.py"
