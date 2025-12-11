module fusion_calling:
    meta_wrapper:
        "v8.0.0/meta/bio/star_arriba"
    config:
        config


use rule * from fusion_calling


use rule star_index from fusion_calling with:
    input:
        fasta=rules.get_genome.output,
        gtf=rules.get_annotation.output,
    threads: 12  # suggestion at: https://github.com/alexdobin/STAR/issues/364#issuecomment-760274128
    resources:
        mem_mb=36000,  # suggestion at: https://github.com/alexdobin/STAR/issues/364#issuecomment-760274128


use rule star_align from fusion_calling with:
    input:
        fq1=get_star_reads_input,
        fq2=lambda wc: get_star_reads_input(wc, True),
        idx=rules.star_index.output,
        annotation=rules.get_annotation.output,
    output:
        aln="results/mapped/star/{sample}.bam",
        reads_per_gene="results/mapped/star/{sample}.ReadsPerGene.tsv",
    params:
        # specific parameters to work well with arriba
        extra=lambda wc, input: "--quantMode GeneCounts "
        f"--sjdbGTFfile {input.annotation} "
        f"{get_star_read_group(wc)}"
        "--outSAMtype BAM SortedByCoordinate "
        "--chimSegmentMin 10 "
        "--chimOutType WithinBAM SoftClip "
        "--chimJunctionOverhangMin 10 "
        "--chimScoreMin 1 "
        "--chimScoreDropMax 30 "
        "--chimScoreJunctionNonGTAG 0 "
        "--chimScoreSeparation 1 "
        "--alignSJstitchMismatchNmax 5 "
        "-1 5 5 "
        "--chimSegmentReadGapMax 3 ",
    resources:
        mem_mb=36000,  # suggestion at: https://github.com/alexdobin/STAR/issues/1159#issuecomment-788150448


ARRIBA_FILTERS_TO_TURN_OFF=",".join(
    [
        # read level filters to turn off (https://github.com/suhrig/arriba/wiki/11-Internal-algorithm#read-level-filters)
        "duplicates",
        "homopolymer",
        "mismatches",
        "low_entropy",
        # event-level filters to turn off (https://github.com/suhrig/arriba/wiki/11-Internal-algorithm#event-level-filters)
        "min_support",
        "no_genomic_support",
        "no_coverage",
        "genomic_support",
    ]
)

use rule arriba from fusion_calling with:
    input:
        bam="results/mapped/star/{sample}.bam",
        genome=rules.get_genome.output,
        annotation=rules.get_annotation.output,
    output:
        fusions="results/arriba/{sample}.fusions.tsv",
        discarded="results/arriba/{sample}.fusions.discarded.tsv",
    params:
        genome_build=config["ref"]["build"],
        default_blacklist=True,
        default_known_fusions=True,
        extra="-u " # do not use arriba-internal duplicate marking
            f"-f {ARRIBA_FILTERS_TO_TURN_OFF}",  # turn off the following arriba filters (https://github.com/suhrig/arriba/wiki/11-Internal-algorithm)
    resources:
        mem_mb=lambda wc, input, attempt: input.size_mb * attempt,  # We have usually seen memory usage well below the input.size_mb (which includes the reference data), but also individual samples with peaks beyond it. One retry to double the reserved memory fixed this for all cases we have seen so far.


rule annotate_exons:
    input:
        fusions="results/arriba/{sample}.fusions.tsv",
        annotation=rules.get_annotation.output,
    output:
        "results/arriba/{sample}.fusions.annotated.tsv",
    conda:
        "../envs/arriba.yaml"
    log:
        "logs/annotate_fusions/{sample}.log",
    shell:
        """
        annotate_exon_numbers.sh {input.fusions} {input.annotation} {output} 2> {log}
        """


rule convert_fusions:
    input:
        fasta=rules.get_genome.output,
        fai=genome_fai,
        fusions="results/arriba/{sample}.fusions.annotated.tsv",
    output:
        temp("results/candidate-calls/arriba/{sample}/{sample}.vcf"),
    conda:
        "../envs/arriba.yaml"
    log:
        "logs/convert_fusions/{sample}.log",
    shell:
        """
        convert_fusions_to_vcf.sh {input.fasta} {input.fusions} {output} 2> {log}
        """


rule sort_arriba_calls:
    input:
        "results/candidate-calls/arriba/{sample}/{sample}.vcf",
    output:
        temp("results/candidate-calls/arriba/{sample}/{sample}.bcf"),
    params:
        # Set to True, in case you want uncompressed BCF output
        uncompressed_bcf=False,
        # Extra arguments
        extras="",
    log:
        "logs/bcf-sort/{sample}/{sample}.log",
    resources:
        mem_mb=8000,
    wrapper:
        "v1.21.0/bio/bcftools/sort"


rule bcftools_concat_candidates:
    input:
        calls=get_arriba_group_candidates,
        idx=lambda wc: get_arriba_group_candidates(wc, csi=True),
    output:
        "results/candidate-calls/arriba/{group}/{group}.bcf",
    log:
        "logs/concat_candidates/{group}.log",
    params:
        uncompressed_bcf=False,
        extra="-d exact -a",
    threads: 4
    resources:
        mem_mb=500,
    wrapper:
        "v1.21.0/bio/bcftools/concat"
