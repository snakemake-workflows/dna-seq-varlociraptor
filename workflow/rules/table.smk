rule vembrane_table:
    input:
        bcf="results/merged-calls/{group}.{event}.fdr-controlled.bcf",
    output:
        bcf="results/tables/{group}.{event}.fdr-controlled.tsv"
    conda:
        "../envs/vembrane.yaml"
    params:
        expression=get_vembrane_expression
    log:
        "logs/vembrane-table/{group}.{event}.log"
    shell:
        "vembrane table \"{params.expression}\" {input.bcf} > {output.bcf} 2> {log}" #"1-10**(-INFO['PROB_DENOVO']/10)"


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
