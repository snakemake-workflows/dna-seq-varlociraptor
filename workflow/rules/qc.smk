rule fastqc:
    input:
        get_fastqc_input,
    output:
        html="results/qc/fastqc/{sample}/{unit}.{fq}.html",
        zip="results/qc/fastqc/{sample}/{unit}.{fq}_fastqc.zip",  # the suffix _fastqc.zip is necessary for multiqc to find the file. If not using multiqc, you are free to choose an arbitrary filename
    log:
        "logs/fastqc/{sample}/{unit}.{fq}.log",
    wrapper:
        "v2.10.0/bio/fastqc"


# TODO
rule samtools_idxstats:
    input:
        bam="results/recal/{sample}.dna.bam",
        idx="results/recal/{sample}.dna.bai",
    output:
        "results/qc/{sample}.dna.bam.idxstats",
    log:
        "logs/samtools/idxstats/{sample}.log",
    wrapper:
        "v2.3.2/bio/samtools/idxstats"


# TODO
rule samtools_stats:
    input:
        bam="results/recal/{sample}.{datatype}.bam",
    output:
        "results/qc/{sample}.{datatype}.bam.stats",
    log:
        "logs/samtools/stats/{sample}_{datatype}.log",
    wrapper:
        "v2.3.2/bio/samtools/stats"


rule multiqc:
    input:
        get_fastqc_results,
    output:
        report(
            "results/qc/multiqc/{group}.{datatype}.html",
            category="Quality control",
            caption="../report/multiqc.rst",
            labels={"Sample group": "{group}"},
        ),
    params:
        "--exclude snippy",
    log:
        "logs/multiqc/{group}_{datatype}.log",
    wrapper:
        "v2.10.0/bio/multiqc"
