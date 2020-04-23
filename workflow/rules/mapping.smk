rule map_reads:
    input:
        reads=get_map_reads_input,
        idx = multiext("resources/genome.fasta", ".amb", ".ann", ".bwt", ".pac", ".sa")
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
        bam=get_recalibrate_quality_input,
        bai=lambda w: get_recalibrate_quality_input(w, bai=True),
        ref="resources/genome.fasta",
        ref_dict="resources/genome.dict",
        ref_fai="resources/genome.fasta.fai",
        known="resources/variation.noiupac.vcf.gz",
        tbi="resources/variation.noiupac.vcf.gz.tbi",
    log:
        "logs/gatk/bqsr/{sample}.log"
    output:
        bam=protected("results/recal/{sample}.sorted.bam")
    params:
        extra=config["params"]["gatk"]["BaseRecalibrator"]
    wrapper:
        "0.47.0/bio/gatk/baserecalibrator"


rule map_primers_se:
    input:
        reads=["results/trimmed/adapters/{sample}/{unit}.single.fastq.gz"],
        idx = multiext(config["primers"]["trimming"]["ref"], ".amb", ".ann", ".bwt", ".pac", ".sa")
    output:
        temp("results/mapped/primers/{sample}/{unit}.single.bam")
    params:
        index=lambda w, input: os.path.splitext(input.idx[0])[0],
        extra = "-L 0 -O 1000 -k 10 -T 10"
    threads: 8
    wrapper:
        "0.51.3/bio/bwa/mem"


rule map_primers_pe:
    input:
        reads=expand("results/trimmed/adapters/{{sample}}/{{unit}}.{read}.fastq.gz", read=[1,2]),
        idx = "{}.fai".format(config["primers"]["trimming "]["ref"])
    output:
        temp("results/mapped/primers/{sample}/{unit}.paired.bam")
    params:
        index=lambda w, input: os.path.splitext(input.idx[0])[0],
        extra = "-L 0,0 -U 1000 -O 1000,1000 -k 10 -T 10"
    threads: 8
    wrapper:
        "0.51.3/bio/bwa/mem"


rule samtools_bam2fq_se:
    input:
        "results/filtered/primers/{sample}/{unit}.single.bam"
    output:
        "results/filtered/primers/{sample}/{unit}.single.fq",
    threads: 3
    wrapper:
        "0.51.3/bio/samtools/bam2fq/interleaved"


rule samtools_bam2fq_pe:
    input:
        "results/filtered/primers/{sample}/{unit}.paired.bam"
    output:
        "results/filtered/primers/{sample}/{unit}.1.fq",
        "results/filtered/primers/{sample}/{unit}.2.fq"
    params:
        sort = "-m 4G",
        bam2fq = "-n"
    threads: 3
    wrapper:
        "0.51.3/bio/samtools/bam2fq/separate"