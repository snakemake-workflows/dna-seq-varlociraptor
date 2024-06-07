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


rule map_reads_vg_giraffe:
    input:
        reads=get_map_reads_input,
        idx=pangenome,
    output:
        "results/mapped/vg/{sample}_mapped.bam",
    log:
        "logs/mapped/vg/{sample}.log",
    benchmark:
        "benchmarks/vg_giraffe/{sample}.tsv"
    conda:
        "../envs/vg.yaml"
    threads: 40
    params:
        lambda wc, input: " -f ".join(input.reads),  # potential issue: in case of single end reads, get map_reads_input() returns a string and join() could create a problem.
    shell:
        "vg giraffe -x {input.idx} -f {params} --output-format BAM -t {threads}  > {output} 2> {log}"  #read groups for samples can be added with: --sample {wildcards.sample} --read-group {wildcards.sample}


rule sort_mapped_vg:
    input:
        "results/mapped/vg/{sample}_mapped.bam",
    output:
        "results/mapped/vg/{sample}_sorted.bam",
    log:
        "logs/samtools_sort_vg/{sample}.log",
    threads: 8
    wrapper:
        "v2.3.2/bio/samtools/sort"


# keep only primary chromosomes
# use -f 2 to keep only properly paired reads, the mate of a read that's on a nonprimary chromosome is problematic for AddOrReplaceReadGroups, because we remove
# all the other nonprimary chromosome from the header.


rule keep_only_primary_chr:
    input:
        "results/mapped/vg/{sample}_sorted.bam",
        "results/mapped/vg/{sample}_sorted.bai",
    output:
        bam="results/mapped/vg/{sample}_extracted.bam",
        idx="results/mapped/vg/{sample}_extracted.bai",
    log:
        "logs/samtools_view_primary_chr/{sample}.log",
    benchmark:
        "benchmarks/samtools_view_primary_chr/{sample}.tsv"
    params:
        region="GRCh38.chr1 GRCh38.chr2 GRCh38.chr3 GRCh38.chr4 GRCh38.chr5 GRCh38.chr6 GRCh38.chr7 GRCh38.chr8 GRCh38.chr9 GRCh38.chr10 GRCh38.chr11 GRCh38.chr12 GRCh38.chr13 GRCh38.chr14 GRCh38.chr15 GRCh38.chr16 GRCh38.chr17 GRCh38.chr18 GRCh38.chr19 GRCh38.chr20 GRCh38.chr21 GRCh38.chr22 GRCh38.chrX GRCh38.chrY GRCh38.chrM",
        extra="-f 2",
    threads: 40
    wrapper:
        "v2.0.0/bio/samtools/view"


# modify the header for chromosome names to be compatible with the reference genome that are acquired from ensembl
# first delete all non classical chromosomes including unlocalized, unplaced and EBV chromosomes (delly complains about them being found in the header)
# second remove GRCh38.chr and third convert M to MT (MT in pangenome reference and M in fasta sequence dict)
# the following sed command replaces the first "M" it finds and replaces it with "MT"


rule reheader:
    input:
        "results/mapped/vg/{sample}_extracted.bam",
    output:
        "results/mapped/vg/{sample}_reheadered.bam",
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

# adding read groups is necessary because base recalibration throws errors
# for not being able to find read group information


rule add_rg:
    input:
        "results/mapped/vg/{sample}_reheadered.bam",
    output:
        "results/mapped/vg/{sample}.bam",
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
        "results/mapped/vg/{sample}.bam",
    output:
        "results/mapped/vg/{sample}.bai",
    log:
        "logs/samtools_index_after_vg_addition/{sample}.log",
    benchmark:
        "benchmarks/samtools_index_after_vg_addition/{sample}.tsv"
    threads: 40
    wrapper:
        "v2.3.2/bio/samtools/index"
