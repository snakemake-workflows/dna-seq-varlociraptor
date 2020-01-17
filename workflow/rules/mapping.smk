rule map_reads:
    input:
        reads=get_merged,
        idx=rules.bwa_index.output
    output:
        temp("mapped/{sample}.sorted.bam")
    log:
        "logs/bwa_mem/{sample}.log"
    params:
        index="refs/genome",
        extra=get_read_group,
        sort="samtools",
        sort_order="coordinate"
    threads: 8
    wrapper:
        "0.39.0/bio/bwa/mem"


rule mark_duplicates:
    input:
        "mapped/{sample}.sorted.bam"
    output:
        bam=temp("dedup/{sample}.sorted.bam"),
        metrics="qc/dedup/{sample}.metrics.txt"
    log:
        "logs/picard/dedup/{sample}.log"
    params:
        config["params"]["picard"]["MarkDuplicates"]
    wrapper:
        "0.39.0/bio/picard/markduplicates"


rule recalibrate_base_qualities:
    input:
        bam="dedup/{sample}.sorted.bam",
        bai="dedup/{sample}.sorted.bam.bai",
        ref="refs/genome.fasta",
        ref_dict="refs/genome.dict",
        known="refs/variation.noiupac.vcf.gz",
        tbi="refs/variation.noiupac.vcf.gz.tbi",
    output:
        bam=protected("recal/{sample}.sorted.bam")
    params:
        extra=config["params"]["gatk"]["BaseRecalibrator"]
    log:
        "logs/gatk/bqsr/{sample}.log"
    wrapper:
        "0.47.0/bio/gatk/baserecalibrator"
