rule get_ganon_reference_input:
    input:
        genome=genome,
    output:
        "resources/ganon/reference.input.tsv",
    log:
        "logs/ganon/get_reference_input.log",
    conda:
        "../envs/pystats.yaml"
    params:
        species=species,
    script:
        "../scripts/get_ganon_reference_input.py"


# The selected reference species is always part of the database.
rule ganon_build_reference:
    input:
        genome=genome,
        info="resources/ganon/reference.input.tsv",
    output:
        multiext(
            "resources/ganon/reference",
            ".hibf",
            ".tax",
        ),
    log:
        "logs/ganon/build_reference.log",
    cache: "omit-software"
    conda:
        "../envs/ganon.yaml"
    threads: 8
    params:
        prefix=subpath(output[0], strip_suffix=".hibf"),
    shell:
        "ganon build-custom --input-file {input.info} --input-target file "
        "--taxonomy ncbi --level species --db-prefix {params.prefix} "
        "--threads {threads} --skip-genome-size 2> {log}"


rule ganon_build_groups:
    output:
        multiext(
            "resources/ganon/groups",
            ".hibf",
            ".tax",
        ),
    log:
        "logs/ganon/build_groups.log",
    cache: "omit-software"
    conda:
        "../envs/ganon.yaml"
    threads: 16
    params:
        prefix=subpath(output[0], strip_suffix=".hibf"),
        organisms=get_ganon_organisms(),
    shell:
        "ganon build {params.organisms} --taxonomy ncbi "
        "--db-prefix {params.prefix} --threads {threads} 2> {log}"


rule ganon_classify:
    input:
        reads=get_sample_fastqs_for_taxonomic_profiling,
        hibf_ref="resources/ganon/reference.hibf",
        hibf_other="resources/ganon/groups.hibf",
        tax_ref="resources/ganon/reference.tax",
        tax_other="resources/ganon/groups.tax",
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
        db_ref=subpath(input.hibf_ref, strip_suffix=".hibf"),
        db_other=subpath(input.hibf_other, strip_suffix=".hibf"),
        reads=get_ganon_reads_param,
    shell:
        "ganon classify --db-prefix {params.db_ref} {params.db_other} --hierarchy-labels reference other "
        "{params.reads} --output-prefix {params.output_prefix} "
        "--threads {threads} 2> {log}"


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
        "--output-file {output} 2> {log}"
