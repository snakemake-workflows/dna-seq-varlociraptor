rule bcf_index:
    input:
        "{prefix}.bcf",
    output:
        "{prefix}.bcf.csi",
    log:
        "logs/bcf-index/{prefix}.log",
    conda:
        "../envs/bcftools.yaml"
    shell:
        "bcftools index {input} 2> {log}"


rule bcf_to_vcf_gz:
    input:
        "{prefix}.bcf",
    output:
        "{prefix}.vcf.gz",
    log:
        "logs/bcf-to-vcf/{prefix}.log",
    wrapper:
        "v1.10.0/bio/bcftools/view"


rule bam_index:
    input:
        "{prefix}.bam",
    output:
        "{prefix}.bai",
    log:
        "logs/bam-index/{prefix}.log",
    wrapper:
        "v1.10.0/bio/samtools/index"


rule tabix_known_variants:
    input:
        "resources/{prefix}.{format}.gz",
    output:
        "resources/{prefix}.{format}.gz.tbi",
    log:
        "logs/tabix/{prefix}.{format}.log",
    params:
        get_tabix_params,
    cache: "omit-software"
    wrapper:
        "v1.12.0/bio/tabix/index"
