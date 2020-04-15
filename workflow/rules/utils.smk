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
        "{prefix}.sorted.bam"
    output:
        "{prefix}.sorted.bam.bai"
    log:
        "logs/bam-index/{prefix}.log"
    wrapper:
        "0.39.0/bio/samtools/index"


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
        "0.45.1/bio/tabix"

rule get_covered_regions:
    input:
        ref="resources/genome.fasta",
        bam="results/recal/{sample}.sorted.bam"
        bai="results/recal/{sample}.sorted.bam.bai"
    output:
        temp("results/regions/temp/{sample}.quantized.bed.gz")
    params:
        prefix={output}.strip(".quantized.bed.gz")
    shadow: "shallow"
    log:
        "logs/bam-regions/{sample}.log"
    conda:
        "../envs/mosdepth.yaml"
    shell:
        "mosdepth {params.prefix} {input.bam} -q 1: 2> {log}"

rule unzip_quantized_regions:
    input:
        "results/regions/temp/{sample}.quantized.bed.gz
    output:
        "results/regions/{sample}.bed"
    shell:
        "gzip -d {input} -c > {output}"