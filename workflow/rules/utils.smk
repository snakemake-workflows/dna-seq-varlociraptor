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


rule yara_index:
    input:
        "resources/genome.fasta"
    output:
        "resources/genome.txt.size",
        "resources/genome.txt.limits",
        "resources/genome.txt.concat",
        "resources/genome.rid.concat",
        "resources/genome.rid.limits",
        "resources/genome.sa.len",
        "resources/genome.sa.val",
        "resources/genome.sa.ind",
        "resources/genome.lf.drp",
        "resources/genome.lf.drs",
        "resources/genome.lf.drv",
        "resources/genome.lf.pst"
    log:
        "logs/yara/index.log"
    conda:
        "../envs/yara.yaml"
    shell:
        "yara_indexer {input} &> {log}"


# ToDO Return all best matches. Requires replacing bamToBed-rules
rule map_primers:
    threads: 12
    input:
        reads=[config["primers"]["trimming"]["primers_fa1"], config["primers"]["trimming"]["primers_fa2"]],
        ref="resources/genome.fasta",
        idx=rules.yara_index.output
    output:
        "results/primers/primers.bam"
    params:
        library_error = config["primers"]["trimming"]["library_error"],
        library_len = config["primers"]["trimming"]["library_length"],
        ref_prefix=lambda w, input: os.path.splitext(input.ref)[0]
    log:
        "logs/yara/map_primers.log"
    conda:
        "../envs/yara.yaml"
    shell:
        "yara_mapper -t {threads} -ll {params.library_len} -ld {params.library_error} -o {output} {params.ref_prefix} {input.reads} > {log}"


rule filter_unmapped_primers:
    input:
        "results/primers/primers.bam"
    output:
        "results/primers/primers.filtered.bam"
    params:
        "-b -f 2"
    log:
        "logs/primers/filtered.log"
    wrapper:
        "0.61.0/bio/samtools/view"


rule primer_to_bedpe:
    input:
        "results/primers/primers.filtered.bam"
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
        "results/primers/primers.filtered.bam"
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
