rule get_target_regions:
    """
    Targets can be defined in one or multiple BED files. In the
    case of multiple BED files, these need to be merged.
    """
    input:
        config.get("target_regions", []),
    output:
        "resources/target_regions/target_regions.bed",
    log:
        "logs/regions/target_regions.log",
    conda:
        "../envs/awk_bedtools.yaml"
    shell:
        """
        (cat {input} | sort -k1,1 -k2,2n - | mergeBed -i - | awk \'{{sub("^chr","", $0); print}}\' > {output} \
        && if [[ ! -s {output} ]]; then >&2 echo 'Empty output: target file appears to be invalid'; exit 1; fi) 2> {log}
        """


rule build_sample_regions:
    input:
        bam="results/recal/{sample}.bam",
        bai="results/recal/{sample}.bai",
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
        "v1.12.0/bio/mosdepth"


rule merge_expanded_group_regions:
    input:
        lambda wc: expand(
            "results/regions/{{group}}/{sample}.quantized.bed.gz",
            sample=get_group_samples(wc.group),
        ),
    output:
        "results/regions/{group}.expanded_regions.bed",
    log:
        "logs/regions/{group}_expanded_regions.log",
    conda:
        "../envs/bedtools.yaml"
    shell:
        "zcat {input} | sort -k1,1 -k2,2n - | mergeBed -i - -d 15000 > {output} 2> {log}"


rule merge_covered_group_regions:
    input:
        lambda wc: expand(
            "results/regions/{{group}}/{sample}.quantized.bed.gz",
            sample=get_group_samples(wc.group),
        ),
    output:
        "results/regions/{group}.covered_regions.bed",
    log:
        "logs/regions/{group}_covered_regions.log",
    conda:
        "../envs/bedtools.yaml"
    shell:
        "zcat {input} | sort -k1,1 -k2,2n - | mergeBed -i - > {output} 2> {log}"


rule filter_group_regions:
    input:
        regions="results/regions/{group}.{regions_type}_regions.bed",
        predefined="resources/target_regions/target_regions.bed"
        if "target_regions" in config
        else [],
        fai=genome_fai,
    output:
        "results/regions/{group}.{regions_type}_regions.filtered.bed",
    conda:
        "../envs/awk_bedtools.yaml"
    params:
        chroms=config["ref"]["n_chromosomes"],
        filter_targets=get_filter_targets,
    log:
        "logs/regions/{group}.{regions_type}_regions.filtered.log",
    shell:
        "cat {input.regions} | grep -f <(head -n {params.chroms} {input.fai} | "
        'awk \'{{print "^"$1"\\t"}}\') {params.filter_targets} | sort -k1,1 -k2,2n '
        "> {output} 2> {log}"


rule build_excluded_group_regions:
    input:
        target_regions="results/regions/{group}.expanded_regions.filtered.bed",
        genome_index=genome_fai,
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
