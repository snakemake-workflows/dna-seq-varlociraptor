rule fastqc:
    input:
        get_fastqc_input,
    output:
        html="results/qc/fastqc/{sample}/{unit}.{fq}.html",
        zip="results/qc/fastqc/{sample}/{unit}.{fq}_fastqc.zip",  # the suffix _fastqc.zip is necessary for multiqc to find the file. If not using multiqc, you are free to choose an arbitrary filename
    log:
        "logs/fastqc/{sample}/{unit}.{fq}.log",
    resources:
        mem_mb=1024,
    wrapper:
        "v2.10.0/bio/fastqc"


rule samtools_idxstats:
    input:
        bam=get_sample_bam,
        idx=subpath(get_sample_bam, strip_suffix=".bam", with_suffix=".bai"),
    output:
        "results/qc/{sample}.bam.idxstats",
    log:
        "logs/samtools/idxstats/{sample}.log",
    wrapper:
        "v9.15.0/bio/samtools/idxstats"


rule samtools_stats:
    input:
        bam=get_sample_bam,
    output:
        "results/qc/{sample}.bam.stats",
    log:
        "logs/samtools/stats/{sample}.log",
    wrapper:
        "v9.15.0/bio/samtools/stats"


rule multiqc:
    input:
        get_multiqc_input,
    output:
        report(
            "results/qc/multiqc/{group}.html",
            category="Quality control",
            caption="../report/multiqc.rst",
            labels={"Sample group": "{group}"},
        ),
    log:
        "logs/multiqc/{group}.log",
    params:
        extra="--exclude snippy",
        use_input_files_only=True,
    wrapper:
        "v9.9.0/bio/multiqc"
