rule map_reads_bwa:
    input:
        reads=get_map_reads_input,
        idx=access.random(rules.bwa_index.output),
    output:
        temp("results/mapped/bwa/{sample}.bam"),
    log:
        "logs/bwa_mem/{sample}.log",
    params:
        extra=get_read_group("-R "),
    threads: 8
    wrapper:
        "v3.8.0/bio/bwa/mem"


rule count_sample_kmers:
    input:
        reads=get_map_reads_input,
    output:
        "results/kmers/{sample}.kff",
    params:
        out_file=lambda wc, output: os.path.splitext(output[0])[0],
        out_dir=lambda wc, output: os.path.dirname(output[0]),
        mem=lambda wc, resources: resources.mem[:-2],
    conda:
        "../envs/kmc.yaml"
    shadow:
        "minimal"
    log:
        "logs/kmers/{sample}.log",
    threads: max(workflow.cores, 1)
    resources:
        mem="64GB",
    shell:
        "kmc -k29 -m{params.mem} -sm -okff -t{threads} -v @<(ls {input.reads}) {params.out_file} {params.out_dir} &> {log}"


rule create_reference_paths:
    output:
        "resources/reference_paths.txt",
    params:
        build=config["ref"]["build"],
    log:
        "logs/reference/paths.log",
    shell:
        'for chrom in {{1..22}} X Y M; do echo "{params.build}#0#chr$chrom"; done > {output} 2> {log}'


rule map_reads_vg:
    input:
        reads=get_map_reads_input,
        graph=access.random(f"{pangenome_prefix}.gbz"),
        kmers=access.random("results/kmers/{sample}.kff"),
        hapl=access.random(f"{pangenome_prefix}.hapl"),
        paths=access.random("resources/reference_paths.txt"),
    output:
        bam=temp("results/mapped/vg/{sample}.raw.bam"),
        indexes=temp(
            multiext(
                f"{pangenome_prefix}.{{sample}}",
                ".gbz",
                ".dist",
                ".shortread.withzip.min",
                ".shortread.zipcodes",
            )
        ),
    log:
        "logs/mapped/vg/{sample}.log",
    benchmark:
        "benchmarks/vg_giraffe/{sample}.tsv"
    params:
        extra=lambda wc, input: f"--ref-paths {input.paths}",
        sorting="none",
    threads: 64
    wrapper:
        "v6.1.0/bio/vg/giraffe"


rule reheader_mapped_reads:
    input:
        "results/mapped/vg/{sample}.raw.bam",
    output:
        temp("results/mapped/vg/{sample}.reheadered.bam"),
    params:
        build=config["ref"]["build"],
    conda:
        "../envs/samtools.yaml"
    log:
        "logs/reheader/{sample}.log",
    shell:
        "(samtools view {input} -H |"
        " sed -E 's/(SN:{params.build}#0#chr)/SN:/; s/SN:M/SN:MT/' | "
        " samtools reheader - {input} > {output}) 2> {log}"


rule fix_mate:
    input:
        "results/mapped/vg/{sample}.reheadered.bam",
    output:
        temp("results/mapped/vg/{sample}.mate_fixed.bam"),
    log:
        "logs/samtools/fix_mate/{sample}.log",
    threads: 8
    params:
        extra="",
    wrapper:
        "v4.7.2/bio/samtools/fixmate"


# adding read groups is exclusive to vg mapped reads and
# necessary because base recalibration throws errors
# for not being able to find read group information
rule add_read_group:
    input:
        lambda wc: (
            "results/mapped/vg/{sample}.mate_fixed.bam"
            if sample_has_primers(wc)
            else "results/mapped/vg/{sample}.reheadered.bam"
        ),
    output:
        temp("results/mapped/vg/{sample}.bam"),
    log:
        "logs/samtools/add_rg/{sample}.log",
    params:
        read_group=get_read_group(""),
        compression_threads=lambda wildcards, threads: (
            f"-@{threads}" if threads > 1 else ""
        ),
    conda:
        "../envs/samtools.yaml"
    threads: 4
    shell:
        "samtools addreplacerg {input} -o {output} -r {params.read_group} "
        "-w {params.compression_threads} 2> {log}"


rule sort_alignments:
    input:
        "results/mapped/{aligner}/{sample}.bam",
    output:
        temp("results/mapped/{aligner}/{sample}.sorted.bam"),
    log:
        "logs/sort/{aligner}/{sample}.log",
    params:
        extra="",
    threads: 16
    resources:
        mem="8GB",
    wrapper:
        "v5.10.0/bio/samtools/sort"


rule annotate_umis:
    input:
        bam="results/mapped/{aligner}/{sample}.sorted.bam",
        idx="results/mapped/{aligner}/{sample}.sorted.bai",
    output:
        temp("results/mapped/{aligner}/{sample}.annotated.bam"),
    conda:
        "../envs/umi_tools.yaml"
    log:
        "logs/annotate_bam/{aligner}/{sample}.log",
    shell:
        "umi_tools group -I {input.bam} --paired --umi-separator : --output-bam -S {output} &> {log}"


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
        idx=access.random(rules.bwa_index.output),
    output:
        temp("results/consensus/{sample}.consensus.{read_type}.mapped.bam"),
    params:
        index=lambda w, input: os.path.splitext(input.idx[0])[0],
        extra=lambda w: "-C " + get_read_group("-R")(w),
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
    threads: 16
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
