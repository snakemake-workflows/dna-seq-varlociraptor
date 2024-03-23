rule get_sra:
    output:
        "sra/{accession}_1.fastq.gz",
        "sra/{accession}_2.fastq.gz",
    log:
        "logs/get-sra/{accession}.log",
    wrapper:
        "v2.3.2/bio/sra-tools/fasterq-dump"


rule cutadapt_pipe:
    input:
        get_cutadapt_pipe_input,
    output:
        pipe("pipe/cutadapt/{sample}/{unit}.{fq}.{ext}"),
    log:
        "logs/pipe-fastqs/catadapt/{sample}-{unit}.{fq}.{ext}.log",
    wildcard_constraints:
        ext=r"fastq|fastq\.gz",
    threads: 0  # this does not need CPU
    shell:
        "cat {input} > {output} 2> {log}"


rule cutadapt_pe:
    input:
        get_cutadapt_input,
    output:
        fastq1=temp("results/trimmed/{sample}/{unit}_R1.fastq.gz"),
        fastq2=temp("results/trimmed/{sample}/{unit}_R2.fastq.gz"),
        qc="results/trimmed/{sample}/{unit}.paired.qc.txt",
    log:
        "logs/cutadapt/{sample}-{unit}.log",
    params:
        extra=config["params"]["cutadapt"],
        adapters=get_cutadapt_adapters,
    threads: 8
    wrapper:
        "v3.5.3/bio/cutadapt/pe"


rule cutadapt_se:
    input:
        get_cutadapt_input,
    output:
        fastq=temp("results/trimmed/{sample}/{unit}.single.fastq.gz"),
        qc="results/trimmed/{sample}/{unit}.single.qc.txt",
    log:
        "logs/cutadapt/{sample}-{unit}.se.log",
    params:
        extra=config["params"]["cutadapt"],
        adapters=get_cutadapt_adapters,
    threads: 8
    wrapper:
        "v3.5.3/bio/cutadapt/se"


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
