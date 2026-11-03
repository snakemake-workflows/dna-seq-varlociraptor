rule ganon_build_groups:
    output:
        multiext(
            "resources/ganon/groups",
            ".hibf",
            ".tax",
        ),
    log:
        "logs/ganon/build_groups.log",
    cache: True
    conda:
        "../envs/ganon.yaml"
    threads: 64
    params:
        prefix=subpath(output[0], strip_suffix=".hibf"),
        ids=ganon_organism_groups,
        select_arg="--organism-group",
        mode=get_ganon_build_mode(),
    shell:
        "(ganon build {params.select_arg} {params.ids} --taxonomy ncbi {params.mode}"
        " --db-prefix {params.prefix} --threads {threads} && rm -rf {params.prefix}_files) &> {log}"


use rule ganon_build_groups as ganon_build_taxids with:
    output:
        multiext(
            "resources/ganon/taxids",
            ".hibf",
            ".tax",
        ),
    log:
        "logs/ganon/build_taxids.log",
    params:
        ids=ganon_taxids,
        select_arg="--taxid",


rule ganon_classify:
    input:
        reads=get_sample_fastqs_for_taxonomic_profiling,
        hibf=get_ganon_db(".hibf"),
        tax=get_ganon_db(".tax"),
    output:
        rep="results/taxonomic_profiling/{sample}.rep",
        tre="results/taxonomic_profiling/{sample}.tre",
    log:
        "logs/ganon/classify/{sample}.log",
    conda:
        "../envs/ganon.yaml"
    threads: 32
    params:
        output_prefix=subpath(output.rep, strip_suffix=".rep"),
        db_prefixes=lambda wildcards, input: [
            subpath(hibf, strip_suffix=".hibf") for hibf in input.hibf
        ],
        reads=get_ganon_reads_param,
    shell:
        "ganon classify --db-prefix {params.db_prefixes} "
        "{params.reads} --output-prefix {params.output_prefix} "
        "--threads {threads} &> {log}"


rule ganon_table:
    input:
        expand(
            "results/taxonomic_profiling/{sample}.tre",
            sample=samples["sample_name"].unique(),
        ),
    output:
        "results/taxonomic_profiling/all.species.tsv",
    log:
        "logs/ganon/table.log",
    conda:
        "../envs/ganon.yaml"
    shell:
        "ganon table --input {input} --rank species --output-value counts "
        "--output-file {output} &> {log}"
