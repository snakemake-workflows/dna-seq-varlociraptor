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
        single_primer=get_single_primer_flag,
    conda:
        "../envs/fgbio.yaml"
    log:
        "logs/trimming/{sample}.log",
    shell:
        "fgbio TrimPrimers -H -i {input.bam} -p {input.primers} -s {params.sort_order} {params.single_primer} -o {output.trimmed} &> {log}"


rule bowtie_build:
    input:
        genome,
    output:
        directory("resources/bowtie_build/"),
    params:
        prefix=f"resources/bowtie_build/{genome_name}",
    log:
        "logs/bowtie/build.log",
    conda:
        "../envs/bowtie.yaml"
    shell:
        "mkdir {output} & "
        "bowtie-build {input} {params.prefix} &> {log}"


rule bowtie_map:
    input:
        primers=lambda w: get_panel_primer_input(w.panel),
        idx="resources/bowtie_build",
    output:
        "results/primers/{panel}_primers.bam",
    params:
        primers=lambda wc, input: format_bowtie_primers(wc, input.primers),
        prefix=f"resources/bowtie_build/{genome_name}",
        insertsize=get_bowtie_insertsize(),
        primer_format=lambda wc, input: "-f" if input_is_fasta(input.primers) else "",
    log:
        "logs/bowtie/{panel}_map.log",
    conda:
        "../envs/bowtie.yaml"
    shell:
        "bowtie {params.primers} -x {params.prefix} {params.insertsize} {params.primer_format} -S | samtools view -b - > {output} 2> {log}"


rule filter_unmapped_primers:
    input:
        "results/primers/{panel}_primers.bam",
    output:
        "results/primers/{panel}_primers.filtered.bam",
    params:
        extra=get_filter_params,
    log:
        "logs/primers/{panel}_primers_filtered.log",
    wrapper:
        "v1.10.0/bio/samtools/view"


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
