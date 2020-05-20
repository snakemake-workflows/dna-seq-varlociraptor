rule bcf_index:
    input:
        "{prefix}.bcf",
    output:
        "{prefix}.bcf.csi"
    log:
        "logs/bcf-index/{prefix}.log"
    conda:
        "../envs/bcftools.yaml"
    shell:
        "bcftools index {input} 2> {log}"


rule bam_index:
    input:
        "{prefix}.sorted.bam"
    output:
        "{prefix}.sorted.bam.bai"
    log:
        "logs/bam-index/{prefix}.log"
    wrapper:
        "0.56.0/bio/samtools/index"


rule tabix_known_variants:
    input:
        "resources/{prefix}.{format}.gz"
    output:
        "resources/{prefix}.{format}.gz.tbi"
    log:
        "logs/tabix/{prefix}.{format}.log"
    params:
        get_tabix_params
    cache: True
    wrapper:
        "0.56.0/bio/tabix"


rule build_genome_bed:
    input:
        "resources/genome.fasta.fai"
    output:
        "results/regions/genome_regions.bed"
    params:
        chroms=config["ref"]["n_chromosomes"]
    log:
        "logs/regions/genome_regions.log"
    script:
        "../scripts/fasta_generate_genome.py"


rule build_excluded_regions:
    input:
        target_regions=config["calling"]["regions"]["bed"],
        genome_regions="results/regions/genome_regions.bed"
    output:
        "results/regions/excluded_regions.bed"
    log:
         "logs/regions/excluded_regions.log"
    conda:
        "../envs/bedtools.yaml"
    shell:
        "bedtools complement -i {input.target_regions} -g {input.genome_regions} > {output} 2> {log}"
