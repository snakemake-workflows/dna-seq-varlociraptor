import glob

rule get_sra:
    output:
        "sra/{accession}_1.fastq",
        "sra/{accession}_2.fastq"
    log:
        "logs/get-sra/{accession}.log"
    wrapper:
        "0.49.0/bio/sra-tools/fasterq-dump"


def get_raw_fastq(wildcards):
    unit = units.loc[wildcards.sample].loc[wildcards.unit]
    if unit["fq1"].endswith("gz"):
        ending = ".gz"
    else:
        ending = ""

    if pd.isna(unit["fq1"]):
        # SRA sample (always paired-end for now)
        accession = unit["sra"]
        return expand("sra/{accession}_{read}.fastq", accession=accession, read=[1, 2])
    if pd.isna(unit["fq2"]):
        # single end local sample
        return f"pipe/cutadapt/{unit.sample_name}-{unit.unit_name}.fq1.fastq{ending}"
    else:
        # paired end local sample
        return expand(f"pipe/cutadapt/{unit.sample_name}-{unit.unit_name}.{{read}}.fastq{ending}", read=["fq1","fq2"])


def cutadapt_pipe_input(wc):
    files = list(sorted(glob.glob(units.loc[wc.sample].loc[wc.unit, wc.fq])))
    #print(wc, units.loc[wc.sample].loc[wc.unit, wc.fq])
    assert(len(files) > 0)
    return files


rule cutadapt_pipe:
    threads:
        0
    input:
        cutadapt_pipe_input
    output:
        pipe("pipe/cutadapt/{sample}-{unit}.{fq}.fast{ending}")
    shell:
        "cat {input} > {output}"


rule cutadapt_pe:
    input:
        get_raw_fastq
    output:
        fastq1="results/trimmed/{sample}-{unit}.1.fastq.gz",
        fastq2="results/trimmed/{sample}-{unit}.2.fastq.gz",
        qc="results/trimmed/{sample}-{unit}.qc.txt"
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
        get_raw_fastq
    output:
        fastq="results/trimmed/{sample}-{unit}.single.fastq.gz",
        qc="results/trimmed/{sample}-{unit}.qc.txt"
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
        lambda w: expand("results/trimmed/{{sample}}-{unit}.{{read}}.fastq.gz", unit=units.loc[w.sample, "unit_name"])
    output:
        "results/merged/{sample}.{read,(single|1|2)}.fastq.gz"
    shell:
        "cat {input} > {output}"
