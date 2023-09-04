rule fastqc:
    input:
        get_fastqc_input,
    output:
        html="results/qc/fastqc/{sample}/{unit}.{fq}.html",
        zip="results/qc/fastqc/{sample}/{unit}.{fq}_fastqc.zip",  # the suffix _fastqc.zip is necessary for multiqc to find the file. If not using multiqc, you are free to choose an arbitrary filename
    log:
        "logs/fastqc/{sample}/{unit}.{fq}.log",
    wrapper:
        "v1.3.2/bio/fastqc"


rule samtools_idxstats:
    input:
        bam="results/recal/{sample}.bam",
        idx="results/recal/{sample}.bai",
    output:
        "results/qc/{sample}.bam.idxstats",
    log:
        "logs/samtools/idxstats/{sample}.log",
    wrapper:
        "v1.10.0/bio/samtools/idxstats"


rule samtools_stats:
    input:
        bam="results/recal/{sample}.bam",
    output:
        "results/qc/{sample}.bam.stats",
    log:
        "logs/samtools/stats/{sample}.log",
    wrapper:
        "v1.10.0/bio/samtools/stats"


rule somalier_preprocess_known_variants:
    input:
        vcf="resources/variation.vcf.gz",
        ref=genome,
        refidx=genome_fai,
    output:
        "resources/variation.somalier.vcf.gz"
    log:
        "logs/somalier/preprocess_known_variants.log"
    conda:
        "../envs/somalier.yaml"
    shadow: "minimal"
    cache: True
    shell:
        "(somalier find-sites <(bcftools annotate -c INFO/AF:=INFO/MAF {input.vcf} | bcftools norm --check-ref x --fasta-ref {input.ref} -);"
        " mv sites.vcf.gz {output}) 2> {log}"


rule somalier_extract:
    input:
        sites="resources/variation.somalier.vcf.gz",
        bam="results/recal/{sample}.bam",
        ref=genome,
        refidx=genome_fai,
    output:
        "results/somalier/{sample}.somalier",
    log:
        "logs/somalier/extract/{sample}.log"
    params:
        outdir=lambda wildcards, output: os.path.dirname(output[0])
    conda:
        "../envs/somalier.yaml"
    shell:
        "somalier extract -d {params.outdir} --sites {input.sites} -f {input.ref} "
        "{input.bam} 2> {log}"


rule somalier_relate:
    input:
        expand("results/somalier/{sample}.somalier", sample=samples["sample_name"])
    output:
        multiext(
            "results/somalier/relatedness/all",
            ".samples.tsv",
            ".pairs.tsv",
            ".groups.tsv",
            ".html",
        )
    params:
        outdir=lambda wildcards, output: os.path.dirname(output[0])
    log:
        "logs/somalier/relate.log"
    conda:
        "../envs/somalier.yaml"
    shell:
        "somalier relate --output-prefix {params.outdir}/all {input} 2> {log}"


rule multiqc:
    input:
        get_multiqc_results,
    output:
        report(
            "results/qc/multiqc/all.html",
            category="Quality control",
            caption="../report/multiqc.rst",
            labels={"Samples": "all"},
        ),
    params:
        extra="--exclude snippy",
        use_input_files_only=True,
    log:
        "logs/multiqc/all.log",
    wrapper:
        "v2.6.0/bio/multiqc"