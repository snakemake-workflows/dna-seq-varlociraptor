rule get_sra:
    output:
        "sra/{accession}_1.fastq.gz",
        "sra/{accession}_2.fastq.gz",
    log:
        "logs/get-sra/{accession}.log",
    wrapper:
        "v5.0.2/bio/sra-tools/fasterq-dump"


rule fastp_pipe:
    input:
        get_fastp_pipe_input,
    output:
        pipe("pipe/fastp/{sample}/{unit}.{fq}.{ext}"),
    log:
        "logs/pipe-fastqs/fastp/{sample}-{unit}.{fq}.{ext}.log",
    wildcard_constraints:
        ext=r"fastq|fastq\.gz",
    threads: 0  # this does not need CPU
    shell:
        "cat {input} > {output} 2> {log}"


rule fastp_se:
    input:
        sample=lambda wc: get_fastp_input(wc),
    output:
        trimmed=temp("results/trimmed/{sample}/{unit}.single.fastq.gz"),
        html="results/trimmed/{sample}/{unit}.se.html",
        json="results/trimmed/{sample}/{unit}.se.json",
    log:
        "logs/fastp/se/{sample}_{unit}.log",
    params:
        adapters=get_fastp_adapters,
        extra=get_fastp_extra,
    threads: 1
    wrapper:
        "v6.2.0/bio/fastp"


rule fastp_pe:
    input:
        sample=lambda wc: get_fastp_input(wc),
    output:
        trimmed=[
            temp("results/trimmed/{sample}/{unit}_R1.fastq.gz"),
            temp("results/trimmed/{sample}/{unit}_R2.fastq.gz"),
        ],
        html="results/trimmed/{sample}/{unit}.html",
        json="results/trimmed/{sample}/{unit}.json",
    log:
        "logs/fastp/pe/{sample}_{unit}.log",
    params:
        adapters=get_fastp_adapters,
        extra=get_fastp_extra,
    threads: 8
    wrapper:
        "v6.2.0/bio/fastp"


# rule cutadapt_pe:
#     input:
#         get_cutadapt_input,
#     output:
#         fastq1=temp("results/trimmed/{sample}/{unit}_R1.fastq.gz"),
#         fastq2=temp("results/trimmed/{sample}/{unit}_R2.fastq.gz"),
#         qc="results/trimmed/{sample}/{unit}.paired.qc.txt",
#     log:
#         "logs/cutadapt/{sample}-{unit}.log",
#     params:
#         extra=config["params"]["cutadapt"],
#         adapters=get_cutadapt_adapters,
#     threads: 8
#     wrapper:
#         "v3.5.3/bio/cutadapt/pe"


# rule cutadapt_se:
#     input:
#         get_cutadapt_input,
#     output:
#         fastq=temp("results/trimmed/{sample}/{unit}.single.fastq.gz"),
#         qc="results/trimmed/{sample}/{unit}.single.qc.txt",
#     log:
#         "logs/cutadapt/{sample}-{unit}.se.log",
#     params:
#         extra=config["params"]["cutadapt"],
#         adapters=get_cutadapt_adapters,
#     threads: 8
#     wrapper:
#         "v3.5.3/bio/cutadapt/se"


rule merge_trimmed_fastqs:
    input:
        get_trimmed_fastqs,
    output:
        "results/merged/{sample}_{read}.fastq.gz",
    log:
        "logs/merge-fastqs/trimmed/{sample}_{read}.log",
    wildcard_constraints:
        read="single|R1|R2",
    shell:
        "cat {input} > {output} 2> {log}"
