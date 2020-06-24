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
    threads: 12
    input:
        reads=[config["primers"]["trimming"]["primers_fa1"], config["primers"]["trimming"]["primers_fa2"]],
        ref="resources/genome.fasta"
    output:
        "results/mapped/primers.bam"
    log:
        "logs/razers3/primers.log"
    params:
        library_error = config["primers"]["trimming"]["library_error"],
        library_len = config["primers"]["trimming"]["library_length"],
        sort_order="1"
    conda:
        "../envs/razers3.yaml"
    shell:
        "razers3 -i 95 -rr 99 -m 100 -tc {threads} -so {params.sort_order} -ll {params.library_len} -le {params.library_error} -o {output} {input.ref} {input.reads}"


rule primer_to_bedpe:
    input:
        "results/mapped/primers.bam"
    output:
        "results/primers/primers.bedpe"
    log:
        "logs/primers/primers_bedpe.log"
    conda:
        "../envs/bedtools.yaml"
    shell:
        "samtools sort -n {input} | bamToBed -i - -bedpe > {output} 2> {log}"


rule primer_to_bed:
    input:
        "results/mapped/primers.bam"
    output:
        "results/primers/primers.bed"
    log:
        "logs/primers/primers_bed.log"
    conda:
        "../envs/bedtools.yaml"
    shell:
        "samtools sort -n {input} | bamToBed -i - > {output} 2> {log}"

rule build_primer_regions:
    input:
        "results/primers/primers.bedpe"
    output:
        "results/primers/primer_regions.tsv"
    log:
        "logs/primers/build_primer_regions.log"
    script:
        "../scripts/build_primer_regions.py"


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
