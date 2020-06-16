rule get_sra:
    output:
        "sra/{accession}_1.fastq",
        "sra/{accession}_2.fastq"
    log:
        "logs/get-sra/{accession}.log"
    wrapper:
        "0.56.0/bio/sra-tools/fasterq-dump"


rule cutadapt_pipe:
    input:
        get_cutadapt_pipe_input
    output:
        pipe('pipe/cutadapt/{sample}/{unit}.{fq}.{ext}')
    log:
        "logs/pipe-fastqs/catadapt/{sample}-{unit}.{fq}.{ext}.log"
    wildcard_constraints:
        ext=r"fastq|fastq\.gz"
    threads: 0 # this does not need CPU
    shell:
        "cat {input} > {output} 2> {log}"


rule cutadapt_pe:
    input:
        get_cutadapt_input
    output:
        fastq1="results/trimmed/{sample}/{unit}_R1.fastq.gz",
        fastq2="results/trimmed/{sample}/{unit}_R2.fastq.gz",
        qc="results/trimmed/{sample}/{unit}.paired.qc.txt"
    log:
        "logs/cutadapt/{sample}-{unit}.log"
    params:
        others = config["params"]["cutadapt"],
        adapters = lambda w: str(units.loc[w.sample].loc[w.unit, "adapters"]),
    threads: 8
    wrapper:
        "0.59.2/bio/cutadapt/pe"

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
        "0.59.2/bio/cutadapt/se"


rule trim_primers:
    input:
        bams="results/mapped/{sample}.sorted.bam",
        primers="results/primers/primer_regions.tsv",
    output:
        "results/mapped/{sample}.trimmed.bam"
    params:
        sort_order="Coordinate"
    conda:
        "../envs/fgbio.yaml"
    log:
        "logs/trimming/{sample}.log"
    shell:
        "fgbio TrimPrimers -H -i {input.bams} -p {input.primers} -s {params.sort_order} -o {output} > {log}"


rule merge_fastqs:
    input:
        get_fastqs
    output:
        "results/merged/{sample}_{read}.fastq.gz"
    log:
        "logs/merge-fastqs/{sample}_{read}.log"
    wildcard_constraints:
        read="single|R1|R2"
    shell:
        "cat {input} > {output} 2> {log}"
