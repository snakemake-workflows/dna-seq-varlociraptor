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
        "{prefix}.bam.bai"
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
