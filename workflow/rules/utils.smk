rule bcf_index:
    input:
        "{prefix}.bcf",
    output:
        "{prefix}.bcf.csi"
    log:
        "logs/bcf-index/{prefix}.log"
    conda:
        "../envs/bcftools.yaml"
    shell:
        "bcftools index {input} 2> {log}"


rule bam_index:
    input:
        "{prefix}.bam"
    output:
        "{prefix}.bai"
    log:
        "logs/bam-index/{prefix}.log"
    wrapper:
        "0.59.2/bio/samtools/index"


rule tabix_known_variants:
    input:
        "resources/{prefix}.{format}.gz"
    output:
        "resources/{prefix}.{format}.gz.tbi"
    log:
        "logs/tabix/{prefix}.{format}.log"
    params:
        get_tabix_params
    cache: True
    wrapper:
        "0.59.2/bio/tabix"


rule testcase:
    input:
        obs=get_group_observations,
        scenario="results/scenarios/{group}.yaml"
    output:
        directory("resources/testcases/{group}.{caller}/{locus}")
    log:
        "logs/varlociraptor/testcase/{group}.{caller}.{locus}.log"
    params:
        obs=lambda w, input: ["{}={}".format(s, f) for s, f in zip(get_group_aliases(w), input.obs)],
        parent=lambda w, output: os.path.dirname(output[0])
    threads: workflow.cores
    conda:
        "../envs/varlociraptor.yaml"
    shell:
        "varlociraptor "
        "call variants --testcase-prefix {output} --testcase-locus {wildcards.locus} "
        "generic --obs {params.obs} "
        "--scenario {input.scenario} 2> {log}"
