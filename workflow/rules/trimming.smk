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
        fastq1="results/trimmed/adapters/{sample}/{unit}.1.fastq.gz",
        fastq2="results/trimmed/adapters/{sample}/{unit}.2.fastq.gz",
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


rule pipe_ptrimmer_pe:
    input:
        "results/trimmed/adapters/{sample}/{unit}.{read}.fastq.gz"
    output:
        pipe("results/trimmed/adapters/{sample}/{unit}.paired_R{read}.fastq.gz")
    log:
        "logs/pipe-fastqs/catadapt/{sample}-{unit}.{read}.log"
    threads: 0 # this does not need CPU
    shell:
        "cat {input} > {output} 2> {log}"

rule ptrimmer_se:
    input:
        "results/trimmed/adapters/{sample}/{unit}.single_R1.fastq.gz"
    output:
        "results/trimmed/primers/{sample}/{unit}.single_trim_R1.fastq.gz",
    log:
        "logs/ptrimmer/{sample}-{unit}.log"
    params:
        primers=config["primers"]["trimming"]["primers"],
        tmp_fq="results/trimmed/adapters/{sample}/{unit}.single_trim_R1.fq",
        folder=lambda wc, output: os.path.dirname(output)
    conda:
        "../envs/ptrimmer.yaml"
    shell:
        "ptrimmer -f {input} -a {params.primers} -o {params.folder}"
        "gzip -c {params.tmp_fq} > {output}"


rule ptrimmer_pe:
    input:
        r1="results/trimmed/adapters/{sample}/{unit}.paired_R1.fastq.gz",
        r2="results/trimmed/adapters/{sample}/{unit}.paired_R2.fastq.gz",
    output:
        fastq1="results/trimmed/primers/{sample}/{unit}.1.fastq.gz",
        fastq2="results/trimmed/primers/{sample}/{unit}.2.fastq.gz",
    log:
        "logs/ptrimmer/{sample}-{unit}.log"
    params:
        primers=config["primers"]["trimming"]["primers"],
        tmp_fq1="results/trimmed/adapters/{sample}/{unit}.paired_trim_R1.fq",
        tmp_fq2="results/trimmed/adapters/{sample}/{unit}.paired_trim_R2.fq",
        folder=lambda wc, output: os.path.dirname(output.fastq1)
    shell:
        "ptrimmer -f {input.r1} -r {input.r2} -a {params.primers} -o {params.folder}"
        "gzip -c {params.tmp_fq1} > {output.fastq1}"
        "gzip -c {params.tmp_fq2} > {output.fastq2}"


rule merge_fastqs:
    input:
        get_trimmed_fastqs
    output:
        "results/merged/{sample}.{read}.fastq.gz"
    log:
        "logs/merge-fastqs/{sample}.{read}.log"
    wildcard_constraints:
        read="single|1|2"
    shell:
        "cat {input} > {output} 2> {log}"
