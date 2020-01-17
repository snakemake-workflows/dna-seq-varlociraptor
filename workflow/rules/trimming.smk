rule cutadapt_pe:
    input:
        lambda w: units.loc[w.sample].loc[w.unit, "fq1"],
        lambda w: units.loc[w.sample].loc[w.unit, "fq2"]
    output:
        fastq1="trimmed/{sample}-{unit}.1.fastq.gz",
        fastq2="trimmed/{sample}-{unit}.2.fastq.gz",
        qc="trimmed/{sample}-{unit}.qc.txt"
    params:
        others = config["params"]["cutadapt"],
        adapters = lambda w: str(units.loc[w.sample].loc[w.unit, "adapters"]),
    log:
        "logs/cutadapt/{sample}-{unit}.log"
    threads: 8
    wrapper:
        "0.42.0/bio/cutadapt/pe"

rule cutadapt_se:
    input:
        lambda w: units.loc[w.sample].loc[w.unit, "fq1"]
    output:
        fastq="trimmed/{sample}-{unit}.fastq.gz",
        qc="trimmed/{sample}-{unit}.qc.txt"
    params:
        others = config["params"]["cutadapt"],
        adapters_r1 = lambda w: str(units.loc[w.sample].loc[w.unit, "adapters"])
    log:
        "logs/cutadapt/{sample}-{unit}.log"
    threads: 8
    wrapper:
        "0.42.0/bio/cutadapt/se"

rule merge_fastqs:
    input:
        lambda w: expand("trimmed/{{sample}}-{unit}.{{read}}.fastq.gz", unit=units.loc[w.sample, "unit_name"])
    output:
        "merged/{sample}.{read,(1|2)}.fastq.gz"
    run:
        if input[0].endswith(".gz"):
            shell("cat {input} > {output}")
        else:
            shell("cat {input} | gzip > {output}")