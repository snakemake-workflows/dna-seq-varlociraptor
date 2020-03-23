rule get_genome:
    output:
        "results/refs/genome.fasta"
    params:
        species=config["ref"]["species"],
        datatype="dna",
        build=config["ref"]["build"],
        release=config["ref"]["release"]
    log:
        "logs/get-genome.log"
    wrapper:
        "0.45.1/bio/reference/ensembl-sequence"


rule genome_faidx:
    input:
        "results/refs/genome.fasta"
    output:
        "results/refs/genome.fasta.fai"
    log:
        "logs/genome-faidx.log"
    wrapper:
        "0.45.1/bio/samtools/faidx"


rule genome_dict:
    input:
        "results/refs/genome.fasta"
    output:
        "results/refs/genome.dict"
    log:
        "logs/picard/create_dict.log"
    wrapper:
        "0.45.1/bio/picard/createsequencedictionary"


rule get_known_variants:
    input:
        # use fai to annotate contig lengths for GATK BQSR
        fai="results/refs/genome.fasta.fai"
    output:
        vcf="results/refs/variation.vcf.gz"
    params:
        species=config["ref"]["species"],
        release=config["ref"]["release"],
        type="all"
    log:
        "logs/get-known-variants.log"
    wrapper:
        "0.45.1/bio/reference/ensembl-variation"


rule remove_iupac_codes:
    input:
        "results/refs/variation.vcf.gz"
    output:
        "results/refs/variation.noiupac.vcf.gz"
    conda:
        "../envs/rbt.yaml"
    log:
        "logs/fix-iupac-alleles.log"
    shell:
        "rbt vcf-fix-iupac-alleles < {input} | bcftools view -Oz > {output}"


rule tabix_known_variants:
    input:
        "results/refs/{prefix}.vcf.gz"
    output:
        "results/refs/{prefix}.vcf.gz.tbi"
    params:
        "-p vcf"
    log:
        "logs/tabix/{prefix}.log"
    wrapper:
        "0.45.1/bio/tabix"


rule bwa_index:
    input:
        "refs/genome.fasta"
    output:
        multiext("results/refs/genome.fasta", ".amb", ".ann", ".bwt", ".pac", ".sa")
    log:
        "logs/bwa_index.log"
    resources:
        mem_mb=369000
    wrapper:
        "0.45.1/bio/bwa/index"
