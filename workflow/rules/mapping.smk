rule map_reads:
    input:
        reads=get_map_reads_input,
        idx=rules.bwa_index.output
    output:
        temp("results/mapped/{sample}.sorted.bam")
    log:
        "logs/bwa_mem/{sample}.log"
    params:
        index=lambda w, input: os.path.splitext(input.idx[0])[0],
        extra=get_read_group,
        sort="samtools",
        sort_order="coordinate"
    threads: 8
    wrapper:
        "0.39.0/bio/bwa/mem"


rule mark_duplicates:
    input:
        "results/mapped/{sample}.sorted.bam"
    output:
        bam=temp("results/dedup/{sample}.sorted.bam"),
        metrics="results/qc/dedup/{sample}.metrics.txt"
    log:
        "logs/picard/dedup/{sample}.log"
    params:
        config["params"]["picard"]["MarkDuplicates"]
    wrapper:
        "0.39.0/bio/picard/markduplicates"


rule recalibrate_base_qualities:
    input:
        bam="results/dedup/{sample}.sorted.bam",
        bai="results/dedup/{sample}.sorted.bam.bai",
        ref="results/refs/genome.fasta",
        ref_dict="results/refs/genome.dict",
        ref_fai="results/refs/genome.fasta.fai",
        known="results/refs/variation.noiupac.vcf.gz",
        tbi="results/refs/variation.noiupac.vcf.gz.tbi",
    log:
        "logs/gatk/bqsr/{sample}.log"
    output:
        bam=protected("results/recal/{sample}.sorted.bam")
    params:
        extra=config["params"]["gatk"]["BaseRecalibrator"]
    wrapper:
        "0.47.0/bio/gatk/baserecalibrator"
