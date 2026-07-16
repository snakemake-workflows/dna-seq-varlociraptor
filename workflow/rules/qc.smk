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
        bam="results/recal/{sample}.bam",
        idx="results/recal/{sample}.bai",
    output:
        "results/qc/{sample}.bam.idxstats",
    log:
        "logs/samtools/idxstats/{sample}.log",
    wrapper:
        "v2.3.2/bio/samtools/idxstats"


rule samtools_stats:
    input:
        bam="results/recal/{sample}.bam",
    output:
        "results/qc/{sample}.bam.stats",
    log:
        "logs/samtools/stats/{sample}.log",
    wrapper:
        "v2.3.2/bio/samtools/stats"


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


rule somalier_find_sites:
    input:
        "resources/variation.noiupac.vcf.gz",
    output:
        "resources/somalier/sites.vcf.gz",
    log:
        "logs/somalier_find_sites.log",
    conda:
        "../envs/somalier.yaml"
    shell:
        "somalier find-sites --min-AF 0.45 --min-AN 2000 --AF-field MAF --AN-field MAC --output-vcf {output} {input} 2> {log}"


rule somalier_extract:
    input:
        bam="results/recal/{sample}.bam",
        sites="resources/somalier/sites.vcf.gz",
        reference=genome,
    output:
        data="results/somalier/data/{sample}.somalier",
    log:
        "logs/somalier_extract/{sample}.log",
    conda:
        "../envs/somalier.yaml"
    params:
        outdir=subpath(output.data, parent=True),
    shell:
        "somalier extract -d {params.outdir} --sites {input.sites} "
        "-f {input.reference} {input.bam} 2> {log}"


rule somalier_relate:
    input:
        collect(
            "results/somalier/data/{sample}.somalier", sample=samples["sample_name"]
        ),
    output:
        samples="results/somalier/all.samples.tsv",
        pairs="results/somalier/all.pairs.tsv",
        groups="results/somalier/all.groups.tsv",
        html="results/somalier/all.html",
    log:
        "logs/somalier_relate.log",
    conda:
        "../envs/somalier.yaml"
    params:
        outdir=subpath(output.samples, parent=True),
    shell:
        "somalier relate -o {params.outdir}/all {input} 2> {log}"
