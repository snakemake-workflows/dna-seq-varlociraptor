rule build_sample_regions:
    input:
        get_varlociraptor_preprocessing_input,
    output:
        temp("results/regions/{group}/{sample}.target_regions.bed"),
    params:
        chroms=config["ref"]["n_chromosomes"],
    log:
        "logs/regions/samples/{group}/{sample}_target_regions.log",
    conda:
        "../envs/bedtools.yaml"
    shell:
        'bamToBed -i {input} | mergeBed -i - | grep -f <(head -n {params.chroms} resources/genome.fasta.fai | awk \'{{print "^"$1"\\t"}}\') > {output} 2> {log}'


rule merge_group_regions:
    input:
        lambda wc: expand(
            "results/regions/{{group}}/{sample}.target_regions.bed",
            sample=get_group_samples(wc.group),
        ),
    output:
        "results/regions/{group}.target_regions.bed",
    log:
        "logs/regions/{group}_target_regions.log",
    conda:
        "../envs/bedtools.yaml"
    shell:
        "cat {input} | sort -k1,1 -k2,2n - | mergeBed -i - > {output} 2> {log}"


rule build_excluded_regions:
    input:
        target_regions="results/regions/{group}.target_regions.bed",
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
