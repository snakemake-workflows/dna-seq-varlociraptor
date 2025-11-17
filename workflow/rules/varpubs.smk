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
        "varpubs deploy-db --db_path {output} --vcf_paths {input.bcf} --email {params.mail} &> {log}"


rule varpubs_summarize_variants:
    input:
        bcf="results/final-calls/{group}.{event}.variants.fdr-controlled.normal-probs.bcf",
        db_path="results/varpubs/{group}.{event}.duckdb",
        cache=get_varpubs_cache(use_before_update=True)
    output:
        summaries="results/varpubs/{group}.{event}.csv",
        cache="results/varpubs/caches/{group}.{event}.duckdb",
    params:
        llm_url=lookup(dpath="varpubs/llm_url", within=config),
        model=lookup(dpath="varpubs/model", within=config),
        api_key=lookup(dpath="varpubs/api_key", within=config),
        cache=lambda wc, input: f"--cache {input.cache}" if input.cache else ""
    conda:
        "../envs/varpubs.yaml"
    log:
        "logs/varpub/summarize/{group}.{event}.log",
    threads: max(workflow.cores, 1)
    shell:
        "varpubs -v summarize-variants --db_path {input.db_path} --vcf_path {input.bcf} --llm_url {params.llm_url} --model {params.model} --api_key '{params.api_key}' --judges 'therapy related' {params.cache} --output {output.summaries} --output_cache {output.cache} &> {log}"


rule merge_caches:
    input:
        cache=get_varpubs_cache(use_before_update=True),
        event_caches=lambda wc: expand(
            "results/varpubs/caches/{group}.{event}.duckdb",
            group=variants_groups,
            event=get_calling_events("variants")
        ),
    output:
        get_varpubs_cache(use_before_update=False),
    log:
        "logs/varpub/merge_cache.log",
    conda:
        "../envs/varpubs.yaml"
    shell:
        "varpubs update-cache --cache {input.event_caches} --output {output} &> {log}"
