rule map_reads:
    input:
        reads=get_map_reads_input,
        idx=rules.bwa_index.output,
    output:
        temp("results/mapped/bwa/{sample}.bam"),
    log:
        "logs/bwa_mem/{sample}.log",
    params:
        extra=get_read_group,
        sorting=get_map_reads_sorting_params,
        sort_order=lambda wc: get_map_reads_sorting_params(wc, ordering=True),
    threads: 8
    wrapper:
        "v3.8.0/bio/bwa/mem"


rule merge_untrimmed_fastqs:
    input:
        get_untrimmed_fastqs,
    output:
        temp("results/untrimmed/{sample}_{read}.fastq.gz"),
    conda:
        "../envs/fgbio.yaml"
    log:
        "logs/merge-fastqs/untrimmed/{sample}_{read}.log",
    wildcard_constraints:
        read="fq1|fq2",
    shell:
        "cat {input} > {output} 2> {log}"


rule sort_untrimmed_fastqs:
    input:
        "results/untrimmed/{sample}_{read}.fastq.gz",
    output:
        temp("results/untrimmed/{sample}_{read}.sorted.fastq.gz"),
    conda:
        "../envs/fgbio.yaml"
    log:
        "logs/fgbio/sort_fastq/{sample}_{read}.log",
    shell:
        "fgbio SortFastq -i {input} -o {output} 2> {log}"


rule annotate_umis:
    input:
        bam="results/mapped/{aligner}/{sample}.bam",
        umi=get_umi_fastq,
    output:
        pipe("pipe/{aligner}/{sample}.annotated.bam"),
    params:
        extra=get_annotate_umis_params,
    log:
        "logs/fgbio/annotate_bam/{aligner}/{sample}.log",
    wrapper:
        "v3.7.0/bio/fgbio/annotatebamwithumis"


rule sort_annotated_reads:
    input:
        "pipe/{aligner}/{sample}.annotated.bam",
    output:
        temp("results/mapped/{aligner}/{sample}.annotated.bam"),
    log:
        "logs/samtools_sort/{aligner}_{sample}.log",
    threads: 8
    wrapper:
        "v3.7.0/bio/samtools/sort"


rule mark_duplicates:
    input:
        bams=get_markduplicates_input,
    output:
        bam=temp("results/dedup/{sample}.bam"),
        metrics="results/qc/dedup/{sample}.metrics.txt",
    log:
        "logs/picard/dedup/{sample}.log",
    params:
        extra=get_markduplicates_extra,
    resources:
        #https://broadinstitute.github.io/picard/faq.html
        mem_mb=3000,
    wrapper:
        "v2.5.0/bio/picard/markduplicates"


rule calc_consensus_reads:
    input:
        get_consensus_input,
    output:
        consensus_r1=temp("results/consensus/fastq/{sample}.1.fq"),
        consensus_r2=temp("results/consensus/fastq/{sample}.2.fq"),
        consensus_se=temp("results/consensus/fastq/{sample}.se.fq"),
        skipped=temp("results/consensus/{sample}.skipped.bam"),
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
        temp("results/consensus/{sample}.consensus.{read_type}.mapped.bam"),
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
        "v2.3.2/bio/bwa/mem"


rule merge_consensus_reads:
    input:
        "results/consensus/{sample}.skipped.bam",
        "results/consensus/{sample}.consensus.se.mapped.bam",
        "results/consensus/{sample}.consensus.pe.mapped.bam",
    output:
        temp("results/consensus/{sample}.merged.bam"),
    log:
        "logs/samtools_merge/{sample}.log",
    threads: 8
    wrapper:
        "v2.3.2/bio/samtools/merge"


rule sort_consensus_reads:
    input:
        "results/consensus/{sample}.merged.bam",
    output:
        temp("results/consensus/{sample}.bam"),
    log:
        "logs/samtools_sort/{sample}.log",
    threads: 8
    wrapper:
        "v2.3.2/bio/samtools/sort"


# TODO Does not use consensus reads
rule splitncigarreads:
    input:
        bam=lambda wc: (
            "results/dedup/{sample}.bam"
            if is_activated("remove_duplicates")
            else "results/mapped/star/{sample}.bam"
        ),
        ref=genome,
    output:
        "results/split/{sample}.bam",
    log:
        "logs/gatk/splitNCIGARreads/{sample}.log",
    params:
        extra="",  # optional
        java_opts="",  # optional
    resources:
        mem_mb=1024,
    wrapper:
        "v3.1.0/bio/gatk/splitncigarreads"


rule recalibrate_base_qualities:
    input:
        bam=get_recalibrate_quality_input,
        bai=lambda w: get_recalibrate_quality_input(w, bai=True),
        ref=genome,
        ref_dict=genome_dict,
        ref_fai=genome_fai,
        known="resources/variation.noiupac.vcf.gz",
        tbi="resources/variation.noiupac.vcf.gz.tbi",
    output:
        recal_table=temp("results/recal/{sample}.grp"),
    params:
        extra=config["params"]["gatk"]["BaseRecalibrator"],
        java_opts="",
    resources:
        mem_mb=1024,
    log:
        "logs/gatk/baserecalibrator/{sample}.log",
    threads: 8
    wrapper:
        "v1.25.0/bio/gatk/baserecalibratorspark"


ruleorder: apply_bqsr > bam_index


rule apply_bqsr:
    input:
        bam=get_recalibrate_quality_input,
        bai=lambda w: get_recalibrate_quality_input(w, bai=True),
        ref=genome,
        ref_dict=genome_dict,
        ref_fai=genome_fai,
        recal_table="results/recal/{sample}.grp",
    output:
        bam=protected("results/recal/{sample}.bam"),
        bai="results/recal/{sample}.bai",
    log:
        "logs/gatk/gatk_applybqsr/{sample}.log",
    params:
        extra=config["params"]["gatk"]["applyBQSR"],  # optional
        java_opts="",  # optional
    wrapper:
        "v2.3.2/bio/gatk/applybqsr"
