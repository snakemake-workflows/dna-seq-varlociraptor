rule gather_calls:
    input:
        calls=gather.calling("results/calls/{{group}}.{{caller}}.{scatteritem}.bcf"),
    output:
        "results/calls/{group}/{caller}.bcf",
    log:
        "logs/gather-calls/{group}/{caller}.log",
    params:
        uncompressed_bcf=False,
        extra="",
    threads: 4
    resources:
        mem_mb=10,
    wrapper:
        "v5.8.2/bio/bcftools/concat"


def get_predictosaurus_build_aux(wildcards):
    aux = []
    for sample in get_group_samples(wildcards.group):
        aux.append(
            "{alias}=results/observations/{group}/{sample}.{caller}.all.bcf".format(
                alias=samples.loc[samples["sample_name"] == sample]["alias"].iloc[0],
                caller=wildcards.caller,
                group=wildcards.group,
                sample=sample,
            )
        )
    return " ".join(aux)


def get_all_group_observations(wildcards):
    return expand(
        "results/observations/{group}/{sample}.{caller}.all.bcf",
        caller=wildcards.caller,
        group=wildcards.group,
        sample=get_group_samples(wildcards.group),
    )


rule predictosaurus_build:
    input:
        calls="results/calls/{group}.{caller}.bcf",
        obs=get_all_group_observations,
    output:
        "results/impact_graphs/{group}.{caller}.graphs.duckdb",
    params:
        obs_aux=get_predictosaurus_build_aux,
    log:
        "logs/predictosaurus/build/{group}_{caller}.log",
    conda:
        "../envs/predictosaurus.yaml"
    shell:
        """
        predictosaurus build -v --calls {input.calls} --observations {params.obs_aux} --output {output} 2> {log}
        """


rule predictosaurus_process:
    input:
        ref=genome,
        gff="resources/annotation.gff3",
        graph="results/impact_graphs/{group}.freebayes.graphs.duckdb",
    output:
        "results/impact_graphs/{group}.scores.duckdb",
    log:
        "logs/predictosaurus/process/{group}.log",
    conda:
        "../envs/predictosaurus.yaml"
    shell:
        """
        predictosaurus process -v --features {input.gff} --reference {input.ref} --graph {input.graph} --output {output} 2> {log}
        """


rule predictosaurus_plot:
    input:
        "results/impact_graphs/{group}.paths.duckdb",
    output:
        "results/impact_graphs/{group}",
    log:
        "logs/predictosaurus/plot/{group}.log",
    conda:
        "../envs/predictosaurus.yaml"
    shell:
        """
        predictosaurus plot --input {input} --output {output} 2> {log}
        """
