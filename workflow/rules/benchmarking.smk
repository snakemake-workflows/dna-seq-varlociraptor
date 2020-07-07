ruleorder: chm_eval_sample > map_reads

rule chm_eval_sample:
    output:
        bam="results/mapped/chm.sorted.bam",
        bai="results/mapped/chm.sorted.bam.bai"
    params:
        # Optionally only grab the first 100 records.
        # This is for testing, remove next line to grab all records.
        first_n=100
    log:
        "logs/benchmarking/chm-eval-sample.log"
    wrapper:
        "0.63.0/bio/benchmark/chm-eval-sample"


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
