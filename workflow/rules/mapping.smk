rule map_reads:
    input:
        reads=get_map_reads_input,
        idx=rules.bwa_index.output,
    output:
        temp("results/mapped/{sample}.sorted.bam"),
    log:
        "logs/bwa_mem/{sample}.log",
    params:
        index=lambda w, input: os.path.splitext(input.idx[0])[0],
        extra=get_read_group,
        sort="samtools",
        sort_order="coordinate",
    threads: 8
    wrapper:
        "0.56.0/bio/bwa/mem"


# TODO: Concider adding read structure as param
# TODO: Set output as temp
rule annotate_umis:
    input:
        bam="results/mapped/{sample}.sorted.bam",
        umi=lambda wc: units.loc[wc.sample]["umis"],
    output:
        "results/mapped/{sample}.annotated.bam",
    params:
        "",
    log:
        "logs/fgbio/annotate_bam/{sample}.log",
    wrapper:
        "v0.80.1/bio/fgbio/annotatebamwithumis"


rule mark_duplicates:
    input:
        lambda wc: "results/mapped/{sample}.sorted.bam"
        if units.loc[wc.sample]["umis"].isnull()
        else "results/mapped/{sample}.annotated.bam",
    output:
        bam=temp("results/dedup/{sample}.sorted.bam"),
        metrics="results/qc/dedup/{sample}.metrics.txt",
    log:
        "logs/picard/dedup/{sample}.log",
    params:
        extra=lambda wc: "{c} {b}".format(
            c=config["params"]["picard"]["MarkDuplicates"],
            b="" if units.loc[wc.sample]["umis"].isnull() else "BARCODE=RX",
        ),
    wrapper:
        "0.80.1/bio/picard/markduplicates"


rule calc_consensus_reads:
    input:
        get_consensus_input,
    output:
        consensus_r1="results/consensus/fastq/{sample}.1.fq",
        consensus_r2="results/consensus/fastq/{sample}.2.fq",
        consensus_se="results/consensus/fastq/{sample}.se.fq",
        skipped="results/consensus/{sample}.skipped.bam",
    log:
        "logs/consensus/{sample}.log",
    conda:
        "../envs/rbt.yaml"
    shell:
        "rbt collapse-reads-to-fragments bam {input} {output} &> {log}"


rule map_consensus_reads:
    input:
        reads=get_processed_consensus_input,
        idx=rules.bwa_index.output,
    output:
        "results/consensus/{sample}.consensus.{read_type}.mapped.bam",
    params:
        index=lambda w, input: os.path.splitext(input.idx[0])[0],
        extra=lambda w: "-C {}".format(get_read_group(w)),
        sort="samtools",
        sort_order="coordinate",
    wildcard_constraints:
        read_type="pe|se",
    log:
        "logs/bwa_mem/{sample}.{read_type}.consensus.log",
    threads: 8
    wrapper:
        "0.67.0/bio/bwa/mem"


rule merge_consensus_reads:
    input:
        "results/consensus/{sample}.skipped.bam",
        "results/consensus/{sample}.consensus.se.mapped.bam",
        "results/consensus/{sample}.consensus.pe.mapped.bam",
    output:
        "results/consensus/{sample}.merged.bam",
    log:
        "logs/samtools_merge/{sample}.log",
    threads: 8
    wrapper:
        "0.67.0/bio/samtools/merge"


rule sort_consensus_reads:
    input:
        "results/consensus/{sample}.merged.bam",
    output:
        "results/consensus/{sample}.sorted.bam",
    log:
        "logs/samtools_sort/{sample}.log",
    threads: 8
    wrapper:
        "0.67.0/bio/samtools/sort"


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
        recal_table=temp("results/recal/{sample}.grp"),
    params:
        extra=config["params"]["gatk"]["BaseRecalibrator"],
        java_opts="",
    log:
        "logs/gatk/baserecalibrator/{sample}.log",
    threads: 8
    wrapper:
        "0.77.0/bio/gatk/baserecalibratorspark"


ruleorder: apply_bqsr > bam_index


rule apply_bqsr:
    input:
        bam=get_recalibrate_quality_input,
        bai=lambda w: get_recalibrate_quality_input(w, bai=True),
        ref="resources/genome.fasta",
        ref_dict="resources/genome.dict",
        ref_fai="resources/genome.fasta.fai",
        recal_table="results/recal/{sample}.grp",
    output:
        bam=protected("results/recal/{sample}.sorted.bam"),
        bai="results/recal/{sample}.sorted.bai",
    log:
        "logs/gatk/gatk_applybqsr/{sample}.log",
    params:
        extra=config["params"]["gatk"]["applyBQSR"],  # optional
        java_opts="",  # optional
    wrapper:
        "0.77.0/bio/gatk/applybqsr"
