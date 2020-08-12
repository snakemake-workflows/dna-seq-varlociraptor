ruleorder: chm_eval_sample > map_reads

rule chm_eval_sample:
    output:
        bam="resources/chm.bam"
    log:
        "logs/benchmarking/chm-eval-sample.log"
    cache: True
    wrapper:
        "master/bio/benchmark/chm-eval-sample"


rule chm_namesort:
    input:
        "resources/chm.bam"
    output:
        pipe("resources/chm.namesorted.bam")
    params:
        "-n -m 4G"
    log:
        "logs/benchmarking/samtools-namesort.log"
    threads: workflow.cores - 1
    wrapper:
        "0.63.0/bio/samtools/sort"


rule chm_to_fastq:
    input:
        "resources/chm.namesorted.bam"
    output:
        fq1="resources/chm.1.fq.gz",
        fq2="resources/chm.2.fq.gz"
    log:
        "logs/benchmarking/samtools-fastq.log"
    conda:
        "../envs/samtools.yaml"
    shell:
        "samtools fastq {input} -1 {output.fq1} -2 {output.fq2} 2> {log}"

rule chm_eval_kit:
    output:
        directory("resources/benchmarking/chm-eval-kit")
    params:
        # Tag and version must match, see https://github.com/lh3/CHM-eval/releases.
        tag="v0.5",
        version="20180222"
    log:
        "logs/benchmarking/chm-eval-kit.log"
    cache: True
    wrapper:
        "0.63.0/bio/benchmark/chm-eval-kit"


rule chm_eval:
    input:
        kit="resources/benchmarking/chm-eval-kit",
        vcf="results/merged-calls/chm.present.fdr-controlled.bcf"
    output:
        summary="benchmarking/chm.summary", # summary statistics
        bed="benchmarking/chm.err.bed.gz" # bed file with errors
    params:
        extra="",
        build="38"
    log:
        "logs/benchmarking/chm-eval.log"
    wrapper:
        "0.63.0/bio/benchmark/chm-eval"
