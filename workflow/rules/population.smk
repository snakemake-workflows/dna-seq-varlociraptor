if is_activated("population/db"):

    rule clean_population_db:
        input:
            get_population_db(use_before_update=True),
        output:
            temp("results/population_db.cleaned.bcf"),
        params:
            remove_groups=lambda wc: ",".join(variants_groups),
        conda:
            "../envs/bcftools.yaml"
        shell:
            "bcftools view --force-samples --samples ^{params.remove_groups} {input} -Ob > {output}"

    rule population_filter_variants:
        input:
            "results/final-calls/{group}.variants.annotated.bcf",
        output:
            temp("results/population/{group}.variants.filtered.bcf"),
        log:
            "logs/population/{group}.filter.log",
        params:
            events=lookup(dpath="population/db/events", within=config),
            fdr=lookup(dpath="population/db/fdr", within=config),
            alias=lookup(dpath="population/db/alias", within=config),
            keep_fields="^INFO/SVLEN,INFO/SVTYPE,INFO/MATEID,INFO/END,INFO/CIPOS,INFO/CIEND,^FORMAT/AF",
        conda:
            "../envs/varlociraptor.yaml"
        shell:
            "varlociraptor filter-calls control-fdr --mode local-smart {input} "
            "--events {params.events} --fdr {params.fdr} | "
            "bcftools view --samples {params.alias} | "
            "bcftools annotate --remove {params.keep_fields} | "
            "bcftools reheader -s <(echo '{wildcards.group}') | "
            "bcftools view -Ou > {output} 2> {log}"

    rule population_db_update:
        input:
            cleaned_db=get_cleaned_population_db(),
            cleaned_db_idx=get_cleaned_population_db(idx=True),
            bcfs=get_population_bcfs(),
            bcfs_idx=get_population_bcfs(idx=True),
        output:
            get_population_db(use_before_update=False),
        log:
            "logs/population/db_export/update_population_db.log",
        conda:
            "../envs/bcftools.yaml"
        shell:
            "bcftools merge -m none {input.cleaned_db} {input.bcfs} -Ob > {output} 2> {log}"
