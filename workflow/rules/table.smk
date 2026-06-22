rule vembrane_table:
    input:
        bcf=get_vembrane_table_input,
        scenario="results/scenarios/{group}.yaml",
    output:
        bcf="results/tables/{group}/{group}.{event}.{calling_type}.fdr-controlled.tsv",
    log:
        "logs/vembrane-table/{group}.{event}.{calling_type}.log",
    conda:
        "../envs/vembrane.yaml"
    params:
        config=lambda wc, input: get_vembrane_config(wc, input),
    shell:
        'vembrane table --wide --header "{params.config[header]}" "{params.config[expr]}" '
        "{input.bcf} > {output.bcf} 2> {log}"


rule tsv_to_excel:
    input:
        tsv="results/{x}.tsv",
    output:
        xlsx="results/{x}.xlsx",
    log:
        "logs/tsv_to_xlsx/{x}.log",
    conda:
        "../envs/excel.yaml"
    script:
        "../scripts/tsv_to_xlsx.py"
