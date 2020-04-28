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
        "logs/pipe-fastqs/catadapt/{sample}-{unit}.{fq}.{ext}"
    wildcard_constraints:
        ext=r"fastq|fastq\.gz"
    threads: 0 # this does not need CPU
    shell:
        "cat {input} > {output} 2> {log}"


rule cutadapt_pe:
    input:
        get_cutadapt_input
    output:
        fastq1="results/trimmed/adapters/{sample}/{unit}_R1.fastq.gz",
        fastq2="results/trimmed/adapters/{sample}/{unit}_R2.fastq.gz",
        qc="results/trimmed/adapters/{sample}/{unit}.paired.qc.txt"
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
        fastq="results/trimmed/adapters/{sample}/{unit}.single.fastq.gz",
        qc="results/trimmed/adapters/{sample}/{unit}.single.qc.txt"
    log:
        "logs/cutadapt/{sample}-{unit}.log"
    params:
        others = config["params"]["cutadapt"],
        adapters_r1 = lambda w: str(units.loc[w.sample].loc[w.unit, "adapters"])
    threads: 8
    wrapper:
        "0.42.0/bio/cutadapt/se"

#TODO Remove rule and set input of ptrimmer_se to unit_R1 (patch pTrimmer)
rule pipe_ptrimmer_se:
    input:
        "results/trimmed/adapters/{sample}/{unit}.single.fastq.gz"
    output:
        pipe("results/trimmed/adapters/{sample}/{unit}.single_R1.fastq.gz")
    log:
        "logs/pipe-fastqs/ptrimmer/{sample}-{unit}.single.log"
    threads: 0 # this does not need CPU
    shell:
        "cat {input} > {output} 2> {log}"


rule ptrimmer_se:
    input:
        "results/trimmed/adapters/{sample}/{unit}.single_R1.fastq.gz"
    output:
        "results/trimmed/primers/{sample}/{unit}.single.fastq.gz",
    log:
        "logs/ptrimmer/{sample}-{unit}.log"
    params:
        primers=config["primers"]["trimming"]["primers"],

    conda:
        "../envs/ptrimmer.yaml"
    shell:
        "ptrimmer -f {input} -a {params.primers} -o ./ > {log}"
        "gzip -c results/trimmed/adapters/{wildcards.sample}/{wildcards.unit}.single_trim_R1.fq > {output}"


rule ptrimmer_pe:
    input:
        r1="results/trimmed/adapters/{sample}/{unit}_R1.fastq.gz",
        r2="results/trimmed/adapters/{sample}/{unit}_R2.fastq.gz",
    output:
        fastq1="results/trimmed/primers/{sample}/{unit}_R1.fastq.gz",
        fastq2="results/trimmed/primers/{sample}/{unit}_R2.fastq.gz",
    log:
        "logs/ptrimmer/{sample}-{unit}.log"
    params:
        primers=config["primers"]["trimming"]["primers"]
    shell:
        "ptrimmer -f {input.r1} -r {input.r2} -a {params.primers} -o ./ > {log}"
        "gzip -c results/trimmed/adapters/{wildcards.sample}/{wildcards.unit}.paired_trim_R1.fq > {output.fastq1}"
        "gzip -c results/trimmed/adapters/{{wildcards.sample}}/{wildcards.unit}.paired_trim_R2.fq > {output.fastq2}"


rule merge_fastqs:
    input:
        get_trimmed_fastqs
    output:
        "results/merged/{sample}_{read}.fastq.gz"
    log:
        "logs/merge-fastqs/{sample}_{read}.log"
    wildcard_constraints:
        read="single|R1|R2"
    shell:
        "cat {input} > {output} 2> {log}"
