rule trim_primers:
    input:
        bams="results/recal/{sample}.filtered.bam",
        primers="results/primers/primer_regions.tsv",
    output:
        "results/trimmed/{sample}.trimmed.bam",
    params:
        sort_order="Coordinate",
    conda:
        "../envs/fgbio.yaml"
    log:
        "logs/trimming/{sample}.log",
    shell:
        "fgbio TrimPrimers -H -i {input.bams} -p {input.primers} -s {params.sort_order} -o {output} &> {log}"


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
        reads=[
            config["primers"]["trimming"]["primers_fa1"],
            config["primers"]["trimming"]["primers_fa2"],
        ],
        ref="resources/genome.fasta",
        idx=rules.yara_index.output,
    output:
        "results/primers/primers.bam",
    params:
        library_len=(
            "-ll {}".format(config["primers"]["trimming"]["library_length"])
            if config["primers"]["trimming"]["library_length"] > 0
            else ""
        ),
        library_error=(
            "-ld {}".format(config["primers"]["trimming"]["library_error"])
            if config["primers"]["trimming"]["library_error"] > 0
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


rule download_liftover_chain:
    output:
        "resources/{origin}To{target}.chain.gz",
    params:
        lambda wc: "https://hgdownload.soe.ucsc.edu/goldenPath/{o}/liftOver/{o}To{t}.over.chain.gz".format(
            o=origin, t=target
        ),
    log:
        "logs/liftover/download_{origin}To{target}.log",
    shell:
        "wget {params} -O {output}"


rule liftover_target_regions:
    input:
        roi=config["target_regions"]["bed"],
        chain="results/primers/{o}To{t}.liftover.bed".format(
            o=config["target_regions"]["liftover"]["origin"],
            t=config["target_regions"]["liftover"]["target"],
        ),
    output:
        "results/primers/target_regions.liftover.bed",
    log:
        "logs/liftover/liftover.log",
    conda:
        "../envs/liftover.yaml"
    shell:
        "liftOver {input.roi} {input.chain} {output} {log} 2> {log}"


rule build_target_regions:
    input:
        "results/primers/primers.bedpe",
    output:
        "results/primers/target_regions.bed",
    log:
        "logs/primers/build_target_regions.log",
    conda:
        "../envs/pandas.yaml"
    script:
        "../scripts/build_target_regions.py"


rule merge_target_regions:
    input:
        "results/primers/target_regions.bed",
    output:
        "results/primers/target_regions.merged.bed",
    log:
        "logs/primers/merge_target_regions.log",
    conda:
        "../envs/bedtools.yaml"
    shell:
        "sort -k1,1 -k2,2n {input} | mergeBed -i - > {output} 2> {log}"


rule build_excluded_regions:
    input:
        target_regions=get_target_regions(),
        genome_index="resources/genome.fasta.fai",
    output:
        "results/primers/excluded_regions.bed",
    params:
        chroms=config["ref"]["n_chromosomes"],
    log:
        "logs/regions/excluded_regions.log",
    conda:
        "../envs/bedtools.yaml"
    shell:
        "(complementBed -i {input.target_regions} -g <(head "
        "-n {params.chroms} {input.genome_index} | cut "
        "-f 1,2 | sort -k1,1 -k 2,2n) > {output}) 2> {log}"


# filter reads that map outside of the expected primer intervals
rule filter_primerless_reads:
    input:
        bam="results/recal/{sample}.sorted.bam",
        regions="results/primers/primers.bed",
    output:
        "results/recal/{sample}.filtered.bam",
    log:
        "logs/primers/{sample}_filter_reads.log",
    conda:
        "../envs/samtools.yaml"
    shell:
        "samtools view -h -b -L {input.regions} {input.bam} > {output} 2> {log}"
