if lookup(dpath="population/db/activate", within=config):

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
            update(lookup(dpath="population/db/path", within=config)),
        log:
            "logs/population/db_export/{group}.log",
        script:
            "scripts/update_population_db.py"
