if is_activated("population/db"):

    rule population_filter_variants:
        input:
            "results/final-calls/{group}.variants.annotated.bcf",
        output:
            "results/population/{group}.variants.filtered.bcf",
        log:
            "logs/population/{group}.filter.log",
        params:
            events=lookup(dpath="population/db/events", within=config),
            fdr=lookup(dpath="population/db/fdr", within=config),
            alias=lookup(dpath="population/db/alias", within=config),
            keep_fields="INFO/SVLEN,INFO/SVTYPE,INFO/MATEID,INFO/END,INFO/CIPOS,INFO/CIEND,FORMAT/AF",
        conda:
            "../envs/varlociraptor.yaml"
        shell:
            "varlociraptor filter-calls control-fdr --mode local-smart {input} "
            "--events {params.events} --fdr {params.fdr} | "
            "bcftools view --samples {params.alias} | "
            "bcftools annotate --remove ^{params.keep_fields} | "
            "bcftools reheader -s <(echo '{wildcards.group}') > {output} 2> {log}"

    rule population_db_update:
        input:
            "results/final-calls/{group}.variants.annotated.bcf",
        output:
            get_population_db(before_update=False),
        log:
            "logs/population/db_export/{group}.log",
        script:
            "scripts/update_population_db.py"

    rule population_db_index:
        input:
            get_population_db(before_update=True),
        output:
            get_population_db(idx=True),
        log:
            "logs/bcf-index/population_db.log",
        conda:
            "../envs/bcftools.yaml"
        shell:
            "bcftools index {input} 2> {log}"
