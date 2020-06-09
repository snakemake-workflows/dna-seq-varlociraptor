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
        "{prefix}.bam"
    output:
        "{prefix}.bam.bai"
    log:
        "logs/bam-index/{prefix}.log"
    wrapper:
        "0.59.2/bio/samtools/index"


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
        "0.59.2/bio/tabix"


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
        sort_order="coordinate",
        extra="-T 10 -k 8 -c 5000"
    threads: 8
    wrapper:
        "0.56.0/bio/bwa/mem"


rule get_primer_insert:
    input:
        "results/mapped/primers.bam",
        "results/mapped/primers.bam.bai"
    output:
        "results/primers/primers.txt"
    log:
        "logs/primers/primers.log"
    conda:
        "../envs/pysam.yaml"
    script:
        "../scripts/extract_primers_insert.py"


rule primer_to_bedpe:
    input:
        "results/mapped/primers.bam"
    output:
        "results/primers/primers.bedpe"
    log:
        "logs/primers/target_regions.log"
    conda:
        "../envs/bedtools.yaml"
    shell:
        "samtools sort -n {input} | bamToBed -i - -bedpe > {output} 2> {log}"

rule build_target_regions:
    input:
        "results/primers/primers.bedpe"
    output:
        "results/primers/target_regions.bed"
    log:
        "logs/primers/build_target_regions.log"
    script:
        "../scripts/build_target_regions.py"


rule merge_target_regions:
    input:
        "results/primers/target_regions.bed"
    output:
        "results/primers/target_regions.merged.bed"
    log:
        "logs/primers/merge_target_regions.log"
    conda:
        "../envs/bedtools.yaml"
    shell:
        "sort -k1,1 -k2,2n {input} | mergeBed -i - > {output} 2> {log}"


rule build_excluded_regions:
    input:
        target_regions="results/primers/target_regions.merged.bed",
        genome_index = "resources/genome.fasta.fai"
    output:
        "results/primers/excluded_regions.bed"
    params:
        chroms=config["ref"]["n_chromosomes"]
    log:
         "logs/regions/excluded_regions.log"
    conda:
        "../envs/bedtools.yaml"
    shell:
        "(complementBed -i {input.target_regions} -g <(head "
        "-n {params.chroms} {input.genome_index} | cut "
        "-f 1,2 | sort -k1,1 -k 2,2n) > {output}) 2> {log}"
