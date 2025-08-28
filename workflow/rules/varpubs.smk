rule varpubs_deploy_db:
    input:
        bcf="results/final-calls/{group}.{event}.variants.fdr-controlled.normal-probs.bcf",
    output:
        "results/varpubs/{group}.{event}.duckdb",
    params:
        mail=lookup(dpath="varpubs/mail", within=config),
    conda:
        "../envs/varpubs.yaml"
    log:
        "logs/varpub/deploy/{group}.{event}.log",
    shell:
        "varpubs deploy-db --db_path {output} --vcf_paths {input.bcf} --email {params.mail}"


rule varpubs_summarize_variants:
    input:
        bcf="results/final-calls/{group}.{event}.variants.fdr-controlled.normal-probs.bcf",
        db_path="results/varpubs/{group}.{event}.duckdb",
    output:
        "results/varpubs/{group}.{event}.csv",
    params:
        llm_url=lookup(dpath="varpubs/llm_url", within=config),
        model=lookup(dpath="varpubs/model", within=config),
        api_key=lookup(dpath="varpubs/api_key", within=config),
    conda:
        "../envs/varpubs.yaml"
    log:
        "logs/varpub/summarize/{group}.{event}.log",
    shell:
        "varpubs summarize-variants --db_path {input.db_path} --vcf_path {input.bcf} --llm_url {params.llm_url} --model {params.model} --api_key {params.api_key} --output {output}"
