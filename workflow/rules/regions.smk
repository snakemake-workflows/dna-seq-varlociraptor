rule build_sample_regions:
    input:
        bam=get_sample_bam,
    output:
        "results/regions/{group}/{sample}.mosdepth.global.dist.txt",
        "results/regions/{group}/{sample}.quantized.bed.gz",
        summary="results/regions/{group}/{sample}.mosdepth.summary.txt",  # this named output is required for prefix parsing
    log:
        "logs/mosdepth/regions/{group}_{sample}.log",
    params:
        extra="--no-per-base",
        quantize="1:",
    wrapper:
        "0.77.0/bio/mosdepth"


rule merge_group_regions:
    input:
        lambda wc: expand(
            "results/regions/{{group}}/{sample}.quantized.bed.gz",
            sample=get_group_samples(wc.group),
        ),
    output:
        "results/regions/{group}.target_regions.bed",
    log:
        "logs/regions/{group}_target_regions.log",
    conda:
        "../envs/bedtools.yaml"
    shell:
        "zcat {input} | sort -k1,1 -k2,2n - | mergeBed -i - -d 15000 > {output} 2> {log}"


rule filter_group_regions:
    input:
        "results/regions/{group}.target_regions.bed",
    output:
        "results/regions/{group}.target_regions.filtered.bed",
    params:
        chroms=config["ref"]["n_chromosomes"],
    log:
        "logs/regions/{group}.target_regions.filtered.log",
    shell:
        'cat {input} | grep -f <(head -n {params.chroms} resources/genome.fasta.fai | awk \'{{print "^"$1"\\t"}}\') > {output} 2> {log}'


rule build_excluded_group_regions:
    input:
        target_regions="results/regions/{group}.target_regions.filtered.bed",
        genome_index="resources/genome.fasta.fai",
    output:
        "results/regions/{group}.excluded_regions.bed",
    params:
        chroms=config["ref"]["n_chromosomes"],
    log:
        "logs/regions/{group}_excluded_regions.log",
    conda:
        "../envs/bedtools.yaml"
    shell:
        "(complementBed -i {input.target_regions} -g <(head "
        "-n {params.chroms} {input.genome_index} | cut "
        "-f 1,2 | sort -k1,1 -k 2,2n) > {output}) 2> {log}"
