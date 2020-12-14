rule vembrane_table:
    input:
        bcf="results/merged-calls/{group}.{event}.fdr-controlled.bcf",
    output:
        bcf="results/tables/{group}.{event}.fdr-controlled.tsv"
    params:
        expression=get_vembrane_expression
    log:
        "logs/vembrane-table/{group}.{event}.log"
    wrapper:
        "master/bio/vembrane/table"


rule tsv_to_excel:
    input:
        tsv="results/{x}.tsv"
    output:
        xlsx="results/{x}.xlsx"
    conda:
        "../envs/excel.yaml"
    log:
        "logs/tsv_to_xlsx/{x}.log"
    script:
        "../scripts/tsv_to_xlsx.py"
