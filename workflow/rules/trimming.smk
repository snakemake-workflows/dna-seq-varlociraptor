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
        "0.59.2/bio/cutadapt/pe"

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
        "0.59.2/bio/cutadapt/se"

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
        fastq="results/trimmed/adapters/{sample}/{unit}.single_R1.fastq.gz",
        primers="results/primers/primers.txt",
    output:
        "results/trimmed/primers/{sample}/{unit}.single.fastq.gz",
    params:
        output_dir=lambda wc, output: os.path.dirname(output)
    log:
        "logs/ptrimmer/{sample}-{unit}.log"
    conda:
        "../envs/ptrimmer.yaml"
    shell:
        "ptrimmer -s single -f {input.fastq} -a {input.primers} -o {params.output_dir} &> {log} && "
        "gzip -c results/trimmed/primers/{wildcards.sample}/{wildcards.unit}.single_trim_R1.fq > {output} && "
        "rm results/trimmed/primers/{wildcards.sample}/{wildcards.unit}.single_trim_R1.fq"


rule ptrimmer_pe:
    input:
        r1="results/trimmed/adapters/{sample}/{unit}_R1.fastq.gz",
        r2="results/trimmed/adapters/{sample}/{unit}_R2.fastq.gz",
        primers="results/primers/primers.txt",
    output:
        r1="results/trimmed/primers/{sample}/{unit}_R1.fastq.gz",
        r2="results/trimmed/primers/{sample}/{unit}_R2.fastq.gz",
    params:
        output_dir=lambda wc, output: os.path.dirname(output.r1)
    log:
        "logs/ptrimmer/{sample}-{unit}.log"
    conda:
        "../envs/ptrimmer.yaml"
    shell:
        "ptrimmer -s pair -f {input.r1} -r {input.r2} -a {input.primers} -o {params.output_dir} &> {log} && "
        "gzip -c -9 results/trimmed/primers/{wildcards.sample}/{wildcards.unit}_trim_R1.fq > {output.r1} && "
        "gzip -c -9 results/trimmed/primers/{wildcards.sample}/{wildcards.unit}_trim_R2.fq > {output.r2} && "
        "rm results/trimmed/primers/{wildcards.sample}/{wildcards.unit}_trim_R1.fq && "
        "rm results/trimmed/primers/{wildcards.sample}/{wildcards.unit}_trim_R2.fq"


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
