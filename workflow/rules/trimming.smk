rule get_sra:
    output:
        "sra/{accession}_1.fastq",
        "sra/{accession}_2.fastq"
    log:
        "logs/get-sra/{accession}.log"
    wrapper:
        "0.49.0/bio/sra-tools/fasterq-dump"


rule cutadapt_pipe:
    input:
        get_cutadapt_pipe_input
    output:
        pipe('pipe/cutadapt/{sample}/{unit}.{fq}.{ext}')
    log:
        "logs/pipe-fastqs/{sample}-{unit}.{fq}.{ext}"
    wildcard_constraints:
        ext=r"fastq|fastq\.gz"
    threads: 0 # this does not need CPU
    shell:
        "cat {input} > {output} 2> {log}"


rule cutadapt_pe:
    input:
        get_cutadapt_input
    output:
        fastq1="results/trimmed/{sample}/{unit}.1.fastq.gz",
        fastq2="results/trimmed/{sample}/{unit}.2.fastq.gz",
        qc="results/trimmed/{sample}/{unit}.paired.qc.txt"
    log:
        "logs/cutadapt/{sample}-{unit}.log"
    params:
        others = config["params"]["cutadapt"],
        adapters = lambda w: str(units.loc[w.sample].loc[w.unit, "adapters"]),
    threads: 8
    wrapper:
        "0.42.0/bio/cutadapt/pe"

rule cutadapt_se:
    input:
        get_cutadapt_input
    output:
        fastq="results/trimmed/{sample}/{unit}.single.fastq.gz",
        qc="results/trimmed/{sample}/{unit}.single.qc.txt"
    log:
        "logs/cutadapt/{sample}-{unit}.log"
    params:
        others = config["params"]["cutadapt"],
        adapters_r1 = lambda w: str(units.loc[w.sample].loc[w.unit, "adapters"])
    threads: 8
    wrapper:
        "0.42.0/bio/cutadapt/se"

rule merge_fastqs:
    input:
        lambda w: expand("results/trimmed/{{sample}}/{unit}.{{read}}.fastq.gz", unit=units.loc[w.sample, "unit_name"])
    output:
        "results/merged/{sample}.{read}.fastq.gz"
    log:
        "logs/merge-fastqs/{sample}.{read}.log"
    wildcard_constraints:
        read="single|1|2"
    shell:
        "cat {input} > {output} 2> {log}"
