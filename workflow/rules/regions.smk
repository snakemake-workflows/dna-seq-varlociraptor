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


rule transform_gene_annotations:
    input:
        "resources/annotation.gtf",
    output:
        "resources/gene_annotation.bed",
    log:
        "logs/plot_coverage/transform_gene_regions.log",
    script:
        "../scripts/transform_gene_regions.py"


rule build_sample_regions:
    input:
        bam="results/recal/{sample}.bam",
        bai="results/recal/{sample}.bai",
        bed="resources/gene_annotation.bed",
    output:
        "results/regions/{group}/{sample}.mosdepth.global.dist.txt",
        "results/regions/{group}/{sample}.quantized.bed.gz",
        "results/regions/{group}/{sample}.regions.bed.gz",
        "results/regions/{group}/{sample}.mosdepth.region.dist.txt",
        summary="results/regions/{group}/{sample}.mosdepth.summary.txt",  # this named output is required for prefix parsing
    log:
        "logs/mosdepth/regions/{group}_{sample}.log",
    params:
        extra="--no-per-base",
        quantize="1:",
    wrapper:
        "v2.3.2/bio/mosdepth"


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
        predefined=(
            "resources/target_regions/target_regions.bed"
            if "target_regions" in config
            else []
        ),
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


rule download_delly_excluded_regions:
    output:
        "results/regions/{species}.{build}.delly_excluded.bed",
    params:
        url="https://raw.githubusercontent.com/dellytools/delly/main/excludeTemplates/{species}.{build}.excl.tsv",
    log:
        "logs/download_delly_regions/{species}_{build}.log",
    shell:
        "curl {params.url} -o {output} &> {log}"
