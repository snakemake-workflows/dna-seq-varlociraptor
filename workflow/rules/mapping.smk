# mapping_mark_duplicates = config["mapping"].get("mark_duplicates", "picard")
# if mapping_mark_duplicates == "picard":
#     mapping_output = "results/mapped/{sample}.cram"
# elif mapping_mark_duplicates == "samblaster":
#     mapping_output = "results/dedup/{sample}.cram"
# else:
#     raise ValueError("Invalid value for duplication tool config: {}".format(mapping_mark_duplicates))

rule map_reads:
    input:
        reads=get_map_reads_input,
        reference=genome,
        idx=rules.bwa_index.output
    output:
         "results/dedup/{sample}.cram"
    log:
        "logs/bwa_meme/{sample}.log",
    params:
        extra=lambda wc: get_read_group(wc) + " -M",
        sort="samtools",  # Can be 'none' or 'samtools'.
        sort_order="coordinate",  # Can be 'coordinate' (default) or 'queryname'.
        sort_extra="",  # Extra args for samtools.
        dedup="mark",  # Can be 'none' (default), 'mark' or 'remove'.
        dedup_extra="-M",  # Extra args for samblaster.
        exceed_thread_limit=True,  # Set threads als for samtools sort / view (total used CPU may exceed threads!)
        embed_ref=True,  # Embed reference when writing cram.
    threads: 88
    wrapper:
        "v1.14.0/bio/bwa-meme/mem"


# ruleorder: map_reads > mark_duplicates

rule annotate_umis:
    input:
        bam="results/dedup/{sample}.cram",
        umi=lambda wc: units.loc[wc.sample]["umis"][0],
    output:
        temp("results/mapped/{sample}.annotated.cram"), # TODO: the reference is missing for the cram 
    resources:
        mem_gb="10",
    log:
        "logs/fgbio/annotate_bam/{sample}.log",
    wrapper:
        "v1.2.0/bio/fgbio/annotatebamwithumis"


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
        "v1.10.0/bio/bwa/mem"


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
        "v1.10.0/bio/samtools/merge"


# rule sort_consensus_reads:
#     input:
#         "results/consensus/{sample}.merged.bam",
#     output:
#         temp("results/consensus/{sample}.bam"),
#     log:
#         "logs/samtools_sort/{sample}.log",
#     threads: 8
#     wrapper:
#         "v1.10.0/bio/samtools/sort"


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
    log:
        "logs/gatk/baserecalibrator/{sample}.log",
    threads: 8
    wrapper:
        "v1.2.0/bio/gatk/baserecalibratorspark"


ruleorder: apply_bqsr > cram_index


rule apply_bqsr:
    input:
        bam=get_recalibrate_quality_input,
        bai=lambda w: get_recalibrate_quality_input(w, bai=True),
        ref=genome,
        ref_dict=genome_dict,
        ref_fai=genome_fai,
        recal_table="results/recal/{sample}.grp",
    output:
        bam="results/recal/{sample}.cram",
    log:
        "logs/gatk/gatk_applybqsr/{sample}.log",
    params:
        extra=config["params"]["gatk"]["applyBQSR"] + " --create-output-bam-index false",  # optional
        spark_master="local[44]",  # optional
        embed_ref=True,
        exceed_thread_limit=True,
    resources:
        mem_mb=10000,
    threads:
        44
    wrapper:
        "v1.16.0/bio/gatk/applybqsrspark"
