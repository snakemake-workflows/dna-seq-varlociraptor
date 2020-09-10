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
        "0.56.0/bio/bwa/mem"


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
        "0.59.2/bio/picard/markduplicates"


rule recalibrate_base_qualities:
    input:
        bam=get_recalibrate_quality_input,
        bai=lambda w: get_recalibrate_quality_input(w, bai=True),
        ref="resources/genome.fasta",
        ref_dict="resources/genome.dict",
        ref_fai="resources/genome.fasta.fai",
        known="resources/variation.noiupac.vcf.gz",
        tbi="resources/variation.noiupac.vcf.gz.tbi",
    output:
        recal_table=temp("results/recal/{sample}.grp")
    params:
        extra=config["params"]["gatk"]["BaseRecalibrator"],
        java_opts=""
    log:
        "logs/gatk/baserecalibrator/{sample}.log"
    threads: 8
    wrapper:
        "0.62.0/bio/gatk/baserecalibratorspark"


ruleorder: apply_bqsr > bam_index


rule apply_bqsr:
    input:
        bam=get_recalibrate_quality_input,
        bai=lambda w: get_recalibrate_quality_input(w, bai=True),
        ref="resources/genome.fasta",
        ref_dict="resources/genome.dict",
        ref_fai="resources/genome.fasta.fai",
        recal_table="results/recal/{sample}.grp"
    output:
        bam=protected("results/recal/{sample}.sorted.bam"),
        bai="results/recal/{sample}.sorted.bai"
    log:
        "logs/gatk/gatk_applybqsr/{sample}.log"
    params:
        extra=config["params"]["gatk"]["applyBQSR"],  # optional
        java_opts="", # optional
    wrapper:
        "0.62.0/bio/gatk/applybqsr"
