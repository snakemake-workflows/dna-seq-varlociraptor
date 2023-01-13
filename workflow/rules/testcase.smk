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
        "v1.16.0/bio/bcftools/concat"


rule testcase:
    input:
        obs=get_all_group_observations,
        scenario="results/scenarios/{group}.yaml",
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
