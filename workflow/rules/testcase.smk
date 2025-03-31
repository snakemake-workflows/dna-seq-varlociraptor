rule sort_observations:
    input:
        "results/observations/{group}/{sample}.{caller}.{scatteritem}.bcf",
    output:
        temp("results/observations/{group}/{sample}.{caller}.{scatteritem}.sorted.bcf"),
    params:
        # Set to True, in case you want uncompressed BCF output
        uncompressed_bcf=False,
        # Extra arguments
        extras="",
    log:
        "logs/bcf-sort/{group}/{sample}.{caller}.{scatteritem}.log",
    resources:
        mem_mb=8000,
    wrapper:
        "v5.5.0/bio/bcftools/sort"


rule index_observations:
    input:
        "results/observations/{group}/{sample}.{caller}.{scatteritem}.sorted.bcf",
    output:
        temp(
            "results/observations/{group}/{sample}.{caller}.{scatteritem}.sorted.bcf.csi"
        ),
    params:
        # Set to True, in case you want uncompressed BCF output
        uncompressed_bcf=False,
        # Extra arguments
        extras="",
    log:
        "logs/bcf-index/{group}/{sample}.{caller}.{scatteritem}.log",
    resources:
        mem_mb=8000,
    wrapper:
        "v5.5.0/bio/bcftools/index"


rule gather_observations:
    input:
        calls=gather.calling(
            "results/observations/{{group}}/{{sample}}.{{caller}}.{scatteritem}.bcf"
        ),
    output:
        "results/observations/{group}/{sample}.{caller}.all.bcf",
    log:
        "logs/gather-observations/{group}/{sample}/{caller}.log",
    params:
        uncompressed_bcf=False,
        extra="",
    threads: 4
    resources:
        mem_mb=10,
    wrapper:
        "v2.3.2/bio/bcftools/concat"


rule testcase:
    input:
        obs=get_all_group_observations,
        scenario="results/scenarios/{group}.yaml",
        ref=genome,
        ref_idx=genome_fai,
        bams=get_group_bams,
        bais=partial(get_group_bams, bai=True),
    output:
        directory("results/testcases/{group}/{caller}/{locus}"),
    log:
        "logs/varlociraptor/testcase/{group}/{caller}/{locus}.log",
    params:
        obs=get_varlociraptor_obs_args,
        extra=config["params"]["varlociraptor"]["call"],
    conda:
        "../envs/varlociraptor.yaml"
    shell:
        "varlociraptor "
        "call variants --testcase-prefix {output} --testcase-locus {wildcards.locus} "
        "{params.extra} generic --obs {params.obs} "
        "--scenario {input.scenario} 2> {log}"
