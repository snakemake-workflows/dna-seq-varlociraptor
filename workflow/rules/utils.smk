rule bcf_index:
    input:
        "{prefix}.bcf",
    output:
        "{prefix}.bcf.csi"
    conda:
        "../envs/bcftools.yaml"
    shell:
        "bcftools index {input}"


rule bam_index:
    input:
        "{prefix}.bam"
    output:
        "{prefix}.bam.bai"
    wrapper:
        "0.39.0/bio/samtools/index"