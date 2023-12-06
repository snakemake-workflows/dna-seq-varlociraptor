rule map_reads:
    input:
        reads=get_map_reads_input,
        idx=rules.bwa_index.output,
    output:
        temp("results/mapped/{sample}.bam"),
    log:
        "logs/bwa_mem/{sample}.log",
    params:
        extra=get_read_group,
        sorting="samtools",
        sort_order="coordinate",
    threads: 8
    wrapper:
        "v2.3.2/bio/bwa/mem"


rule merge_untrimmed_fastqs:
    input:
        get_untrimmed_fastqs,
    output:
        temp("results/untrimmed/{sample}_{read}.fastq.gz"),
    log:
        "logs/merge-fastqs/untrimmed/{sample}_{read}.log",
    wildcard_constraints:
        read="fq1|fq2",
    shell:
        "cat {input} > {output} 2> {log}"


rule annotate_umis:
    input:
        bam=get_mapped_input,
        umi=get_umi_fastq,
    output:
        temp("results/mapped/{sample}.annotated.bam"),
    params:
        extra=get_umi_read_structure,
    resources:
        mem_mb=lambda wc, input: 2.5 * input.size_mb,
    log:
        "logs/fgbio/annotate_bam/{sample}.log",
    wrapper:
        "v2.3.2/bio/fgbio/annotatebamwithumis"


rule mark_duplicates:
    input:
        bams=lambda wc: "results/mapped/{sample}.annotated.bam"
        if sample_has_umis(wc.sample)
        else get_mapped_input(wc),
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


rule map_reads_vg_giraffe:
    input:
        reads=get_map_reads_input,
        idx="resources/pangenome/hprc-v1.0-mc-grch38.xg",
    output:
        "results/vg_mapped/{sample}.bam",
    log:
        "logs/vg_mapped/{sample}.log",
    benchmark:
        "benchmarks/vg_giraffe/{sample}.tsv"
    conda:
        "../envs/vg.yaml"
    threads: 40
    params:
        lambda wc, input: " -f ".join(input.reads),  # potential issue: in case of single end reads, get map_reads_input() returns a string and join() could create a problem.
    shell:
        "vg giraffe -x {input.idx} -f {params} --output-format BAM -t {threads}  > {output} 2> {log}"


rule sort_vg_mapped:
    input:
        "results/vg_mapped/{sample}.bam",
    output:
        "results/vg_mapped/{sample}_sorted.bam",
    log:
        "logs/samtools_sort_vg/{sample}.log",
    threads: 8
    wrapper:
        "v2.3.2/bio/samtools/sort"

#keep only primary chromosomes
rule keep_only_primary_chr:
    input:
        "results/vg_mapped/{sample}_sorted.bam",
        "results/vg_mapped/{sample}_sorted.bai",
    output:
        bam="results/vg_mapped/{sample}_extracted.bam",
        idx="results/vg_mapped/{sample}_extracted.bai"
    log:
        "logs/samtools_view_primary_chr/{sample}.log",
    benchmark:    
        "benchmarks/samtools_view_primary_chr/{sample}.tsv"
    params:
        region="1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 X Y MT"
    threads: 40
    wrapper:
        "v2.0.0/bio/samtools/view"

# modify the header for chromosome names to be compatible with the reference genome that are acquired from ensembl
# first delete all non classical chromosomes including unlocalized, unplaced and EBV chromosomes (delly complains about them being found in the header)
# second remove GRCh38.chr and third convert M to MT (MT in pangenome reference and M in fasta sequence dict)
# the following sed command replaces the first "M" it finds and replaces it with "MT"

rule reheader:
    input:
        "results/vg_mapped/{sample}_extracted.bam",
    output:
        "results/vg_mapped/{sample}_reheadered.bam",
    log:
        "logs/samtools_reheader/{sample}.log",
    benchmark:
        "benchmarks/samtools_reheader/{sample}.tsv"
    conda:
        "../envs/samtools.yaml"
    threads: 40
    shell:
        "samtools view -H {input} | sed '/random/d;/chrUn/d;/EBV/d;s/GRCh38.chr//g;0,/M/s//MT/' | samtools reheader - {input} > {output} 2> {log}"

# rule samtools_index_after_reheader:
#     input:
#         "results/vg_mapped/{sample}_reheadered.bam",
#     output:
#         "results/vg_mapped/{sample}_reheadered.bai",
#     log:
#         "logs/samtools_index_after_reheader/{sample}.log",
#     benchmark:
#         "benchmarks/samtools_index_after_reheader/{sample}.tsv"
#     threads: 40
#     wrapper:
#         "v2.3.2/bio/samtools/index"

#adding read groups is necessary because base recalibration throws errors 
#for not being able to find read group information

rule add_rg:
    input:
        "results/vg_mapped/{sample}_reheadered.bam",
    output:
        "results/vg_mapped/{sample}_rg_added.bam",
    log:
        "logs/picard/add_rg/{sample}.log",
    params:
        extra="--RGLB lib1 --RGPL illumina --RGPU {sample} --RGSM {sample}",
    resources:
        mem_mb=60000,
    wrapper:
        "v2.3.2/bio/picard/addorreplacereadgroups"

rule samtools_index_after_rg_addition:
    input:
        "results/vg_mapped/{sample}_rg_added.bam",
    output:
        "results/vg_mapped/{sample}_rg_added.bai",
    log:
        "logs/samtools_index_after_vg_addition/{sample}.log",
    benchmark:
        "benchmarks/samtools_index_after_vg_addition/{sample}.tsv"
    threads: 40
    wrapper:
        "v2.3.2/bio/samtools/index"