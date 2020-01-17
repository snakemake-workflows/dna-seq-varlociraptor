rule get_genome:
    output:
        "refs/genome.fasta"
    params:
        species=config["ref"]["species"],
        datatype="dna",
        build=config["ref"]["build"],
        release=config["ref"]["release"]
    wrapper:
        "0.45.1/bio/reference/ensembl-sequence"


rule genome_faidx:
    input:
        "refs/genome.fasta"
    output:
        "refs/genome.fasta.fai"
    wrapper:
        "0.45.1/bio/samtools/faidx"


rule genome_dict:
    input:
        "refs/genome.fasta"
    output:
        "refs/genome.dict"
    log:
        "logs/picard/create_dict.log"
    wrapper:
        "0.45.1/bio/picard/createsequencedictionary"


rule get_known_variants:
    input:
        # use fai to annotate contig lengths for GATK BQSR
        fai="refs/genome.fasta.fai"
    output:
        vcf="refs/variation.vcf.gz"
    params:
        species=config["ref"]["species"],
        release=config["ref"]["release"],
        type="all"
    wrapper:
        "0.45.1/bio/reference/ensembl-variation"


rule remove_iupac_codes:
    input:
        "refs/variation.vcf.gz"
    output:
        "refs/variation.noiupac.vcf.gz"
    conda:
        "../envs/rbt.yaml"
    shell:
        "rbt vcf-fix-iupac-alleles < {input} | bcftools view -Oz > {output}"


rule tabix_known_variants:
    input:
        "refs/{prefix}.vcf.gz"
    output:
        "refs/{prefix}.vcf.gz.tbi"
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
        multiext("refs/genome", ".amb", ".ann", ".bwt", ".pac", ".sa")
    log:
        "logs/bwa_index.log"
    params:
        prefix="refs/genome"
    wrapper:
        "0.45.1/bio/bwa/index"
