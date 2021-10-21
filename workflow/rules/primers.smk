rule assign_primers:
    input:
        bam=get_trimming_input,
        primers=get_primer_regions,
    output:
        assigned="results/primers/{sample}.assigned.bam",
        metric="results/primers/{sample}.metric.bam",
    conda:
        "../envs/fgbio.yaml"
    log:
        "logs/primers/assignment/{sample}.log",
    shell:
        "fgbio AssignPrimers -i {input.bam} -p {input.primers} -m {output.metric} -o {output.assigned} &> {log}"


rule filter_primerless_reads:
    input:
        "results/primers/{sample}.assigned.bam",
    output:
        primers="results/primers/{sample}.primers.bam",
        primerless="results/primers/{sample}.primerless.bam",
    conda:
        "../envs/filter_reads.yaml"
    log:
        "logs/primers/filter/{sample}.log",
    script:
        "../scripts/filter_primers.rs"


rule trim_primers:
    input:
        bam="results/primers/{sample}.primers.bam",
        primers=get_primer_regions,
    output:
        trimmed="results/trimmed/{sample}.trimmed.bam",
    params:
        sort_order="Coordinate",
        single_primer=lambda w: (
            "--first-of-pair"
            if not isinstance(get_sample_primer_fastas(w.sample), list)
            else ""
        ),
    conda:
        "../envs/fgbio.yaml"
    log:
        "logs/trimming/{sample}.log",
    shell:
        "fgbio TrimPrimers -H -i {input.bam} -p {input.primers} -s {params.sort_order} {params.single_primer} -o {output.trimmed} &> {log}"


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
        get_primer_bed,
    output:
        "results/primers/{panel}_primer_regions.tsv",
    log:
        "logs/primers/build_{panel}_primer_regions.log",
    conda:
        "../envs/pandas.yaml"
    script:
        "../scripts/build_primer_regions.py"
