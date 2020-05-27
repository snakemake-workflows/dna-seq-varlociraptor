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


rule map_primers:
    input:
        reads=[config["primers"]["trimming"]["primers_fa1"], config["primers"]["trimming"]["primers_fa2"]],
        idx=rules.bwa_index.output
    output:
        "results/mapped/primers.bam"
    log:
        "logs/bwa_mem/primers.log"
    params:
        index=lambda w, input: os.path.splitext(input.idx[0])[0],
        sort="samtools",
        sort_order="queryname",
        extra="-T 10 -k 8 -c 5000"

    threads: 8
    wrapper:
        "0.56.0/bio/bwa/mem"


rule get_primer_insert:
    input:
        "results/mapped/primers.bam"
    output:
        "results/primers/primers.txt"
    log:
        "logs/primers/primers.log"
    conda:
        "../envs/bedtools.yaml"
    shell:
        "bedtools bamtobed -i {input} -bedpe| awk '{{print $1 \"\t\" $5-$3}}' > {output}"


rule get_primer_interval:
    input:
        "results/mapped/primers.bam"
    output:
        "results/primers/target_regions.bed"
    log:
        "logs/primers/target_regions.log"
    conda:
        "../envs/bedtools.yaml"
    shell:
        "bedtools bamtobed -i {input} -bedpe| awk '{{print $1 \"\t\" $2 \"\t\" $6}}' > {output}"


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
        target_regions="results/primers/target_regions.bed",
        genome_regions="results/regions/genome_regions.bed"
    output:
        "results/primers/excluded_regions.bed"
    log:
         "logs/regions/excluded_regions.log"
    conda:
        "../envs/bedtools.yaml"
    shell:
        "bedtools complement -i {input.target_regions} -g {input.genome_regions} > {output} 2> {log}"
