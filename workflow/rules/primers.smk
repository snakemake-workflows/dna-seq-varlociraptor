# TODO Wait for merge/Add own jar file
rule trim_primers:
    input:
        bams="results/recal/{sample}.sorted.bam",
        primers=get_primer_regions,
    output:
        trimmed="results/trimmed/{sample}.trimmed.bam",
        primerless="results/trimmed/{sample}.primerless.bam",
    params:
        sort_order="Coordinate",
        single_primer=lambda w: (
            "--first-of-pair"
            if not isinstance(get_sample_primer_fastas(w.sample), list)
            else ""
        ),
    conda:
        #"../envs/fgbio.yaml"
        "../envs/fgbio_tmp.yaml"
    log:
        "logs/trimming/{sample}.log",
    shell:
        "java -jar ../workflow/scripts/fgbio.jar TrimPrimers -H -i {input.bams} -p {input.primers} -s {params.sort_order} {params.single_primer} -o {output.trimmed} -u {output.primerless} &> {log}"


rule bowtie_build:
    input:
        "resources/genome.fasta",
    output:
        directory("resources/bowtie_build/"),
    params:
        prefix="resources/bowtie_build/genome.fasta",
    log:
        "logs/bowtie/build.log",
    conda:
        "../envs/bowtie.yaml"
    shell:
        "mkdir {output} & "
        "bowtie-build {input} {params.prefix} &> {log}"


rule bowtie_map:
    input:
        reads=lambda w: get_panel_primer_fastas(w.panel),
        idx="resources/bowtie_build",
    output:
        "results/primers/{panel}_primers.bam",
    params:
        reads=(
            lambda wc, input: "-1 {r1} -2 {r2}".format(
                r1=input.reads[0], r2=input.reads[1]
            )
            if isinstance(input.reads, list)
            else "-f {}".format(input.reads)
        ),
        prefix="resources/bowtie_build/genome.fasta",
        insertsize=(
            "-X {}".format(config["primers"]["trimming"].get("library_length"))
            if config["primers"]["trimming"].get("library_length", 0) != 0
            else ""
        ),
    log:
        "logs/bowtie/{panel}_map.log",
    conda:
        "../envs/bowtie.yaml"
    shell:
        "bowtie {params.reads} -x {params.prefix} {params.insertsize} -S | samtools view -b - > {output} 2> {log}"


"""
rule yara_index:
    input:
        "resources/genome.fasta",
    output:
        multiext(
            "resources/genome",
            ".txt.size",
            ".txt.limits",
            ".txt.concat",
            ".rid.limits",
            ".rid.concat",
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
"""


# TODO This might be done by separating bowtie output
rule filter_unmapped_primers:
    input:
        "results/primers/{panel}_primers.bam",
    output:
        "results/primers/{panel}_primers.filtered.bam",
    params:
        extra=(
            lambda wc: "-b -f 2"
            if isinstance(get_panel_primer_fastas(wc.panel), list)
            else "-b -F 4"
        ),
    log:
        "logs/primers/{panel}_primers_filtered.log",
    wrapper:
        "0.61.0/bio/samtools/view"


rule primer_to_bed:
    input:
        "results/primers/{panel}_primers.filtered.bam",
    output:
        "results/primers/{panel}_primers.{ext}",
    wildcard_constraints:
        ext="bedpe|bed",
    params:
        format=lambda wc: "-bedpe" if wc.ext == "bedpe" else "",
    log:
        "logs/primers/{panel}_primers_{ext}.log",
    conda:
        "../envs/bedtools.yaml"
    shell:
        "samtools sort -n {input} | bamToBed -i - {params.format} > {output} 2> {log}"


rule build_primer_regions:
    input:
        lambda w: get_primer_bed(w),
    output:
        "results/primers/{panel}_primer_regions.tsv",
    log:
        "logs/primers/build_{panel}_primer_regions.log",
    conda:
        "../envs/pandas.yaml"
    script:
        "../scripts/build_primer_regions.py"
