# TODO Wait for merge/Add own jar file
rule trim_primers:
    input:
        bams="results/recal/{sample}.sorted.bam",
        primers="results/primers/primer_regions.tsv",
    output:
        trimmed="results/trimmed/{sample}.trimmed.bam",
        primerless="results/trimmed/{sample}.primerless.bam",
    params:
        sort_order="Coordinate",
        single_primer=(
            "--first-of-pair"
            if not config["primers"]["trimming"].get("primers_fa2", "")
            else ""
        ),
    conda:
        "../envs/fgbio.yaml"
    log:
        "logs/trimming/{sample}.log",
    shell:
        "fgbio TrimPrimers -H -i {input.bams} -p {input.primers} -s {params.sort_order} -o {output.trimmed} -u {output.primerless} &> {log}"


rule yara_index:
    input:
        "resources/genome.fasta",
    output:
        multiext(
            "resources/genome.fasta",
            ".txt.size",
            ".txt.limits",
            ".txt.concat",
            ".rid.limits",
            ".sa.len",
            ".sa.val",
            ".sa.ind",
            ".lf.drp",
            ".lf.drs",
            ".lf.drv",
            ".lf.pst",
        ),
    log:
        "logs/yara/index.log",
    conda:
        "../envs/yara.yaml"
    shell:
        "yara_indexer {input} &> {log}"


# ToDO Return all best matches. Requires replacing bamToBed-rules
rule map_primers:
    threads: 12
    input:
        reads=get_primer_fastqs(),
        ref="resources/genome.fasta",
        idx=rules.yara_index.output,
    output:
        "results/primers/primers.bam",
    params:
        library_len=(
            "-ll {}".format(config["primers"]["trimming"]["library_length"])
            if config["primers"]["trimming"].get("library_length", 0) > 0
            else ""
        ),
        library_error=(
            "-ld {}".format(config["primers"]["trimming"]["library_error"])
            if config["primers"]["trimming"].get("library_error", 0) > 0
            else ""
        ),
        ref_prefix=lambda w, input: os.path.splitext(input.ref)[0],
    log:
        "logs/yara/map_primers.log",
    conda:
        "../envs/yara.yaml"
    shell:
        "yara_mapper -t {threads} {params.library_len} {params.library_error} -o {output} {params.ref_prefix} {input.reads} > {log}"


rule filter_unmapped_primers:
    input:
        "results/primers/primers.bam",
    output:
        "results/primers/primers.filtered.bam",
    params:
        "-b -f 2",
    log:
        "logs/primers/filtered.log",
    wrapper:
        "0.61.0/bio/samtools/view"


rule primer_to_bed:
    input:
        "results/primers/primers.filtered.bam",
    output:
        "results/primers/primers.{ext}",
    wildcard_constraints:
        ext="bedpe|bed",
    params:
        format=lambda wc: "-bedpe" if wc.ext == "bedpe" else "",
    log:
        "logs/primers/primers_{ext}.log",
    conda:
        "../envs/bedtools.yaml"
    shell:
        "samtools sort -n {input} | bamToBed -i - {params.format} > {output} 2> {log}"


rule build_primer_regions:
    input:
        get_primer_bed(),
    output:
        "results/primers/primer_regions.tsv",
    log:
        "logs/primers/build_primer_regions.log",
    conda:
        "../envs/pandas.yaml"
    script:
        "../scripts/build_primer_regions.py"


rule build_sample_regions:
    input:
        "results/trimmed/{sample}.trimmed.bam",
    output:
        temp("results/regions/samples/{sample}.target_regions.bed"),
    log:
        "logs/regions/samples/{sample}_target_regions.log",
    conda:
        "../envs/bedtools.yaml"
    shell:
        "bamToBed -i {input} | mergeBed -i - > {output} 2> {log}"


rule merge_group_regions:
    input:
        lambda wc: expand(
            "results/regions/samples/{sample}.target_regions.bed",
            sample=get_group_samples(wc.group),
        ),
    output:
        "results/regions/groups/{group}.target_regions.bed",
    log:
        "logs/regions/groups/{group}_target_regions.log",
    conda:
        "../envs/bedtools.yaml"
    shell:
        "cat {input} | sort -k1,1 -k2,2n - | mergeBed -i - > {output} 2> {log}"


rule build_excluded_regions:
    input:
        target_regions="results/regions/groups/{group}.target_regions.bed",
        genome_index="resources/genome.fasta.fai",
    output:
        "results/regions/groups/{group}.excluded_regions.bed",
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
