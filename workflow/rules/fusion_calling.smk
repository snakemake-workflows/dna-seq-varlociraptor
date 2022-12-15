module fusion_calling:
    meta_wrapper:
        "v1.21.0/meta/bio/star_arriba"
    config:
        config


use rule * from fusion_calling


use rule star_index from fusion_calling with:
    input:
        fasta=rules.get_genome.output,
        annotation=rules.get_annotation.output,


use rule star_align from fusion_calling with:
    input:
        fq1=get_star_reads_input,
        fq2=lambda wc: get_star_reads_input(wc, True),
        idx=rules.star_index.output,
        annotation=rules.get_annotation.output,
    output:
        aln="results/mapped_arriba/{sample}.bam",
        reads_per_gene="results/mapped_arriba/{sample}_ReadsPerGene.out.tab",


use rule arriba from fusion_calling with:
    input:
        bam=get_arriba_group_bam,
        genome=rules.get_genome.output,
        annotation=rules.get_annotation.output,
    output:
        fusions="results/arriba/{group}.fusions.tsv",
        discarded="results/arriba/{group}.fusions.discarded.tsv",
    log:
        "logs/arriba/{group}.log",


rule convert_fusions:
    input:
        fasta=rules.get_genome.output,
        fusions="results/arriba/{group}.fusions.tsv",
    output:
        temp("results/candidate-calls/{group}.arriba.vcf"),
    conda:
        "../envs/arriba.yaml"
    log:
        "logs/convert_fusions/{group}.log",
    shell:
        """
        convert_fusions_to_vcf.sh {input.fasta} {input.fusions} {output} &> {log}
        """


rule bcftools_reheader:
    input:
        vcf="results/candidate-calls/{group}.arriba.vcf",
        fai=genome_fai,
    output:
        "results/candidate-calls/{group}.arriba.bcf",
    conda:
        "../envs/bcftools.yaml"
    log:
        "logs/reheader/{group}.log",
    shell:
        """
        bcftools reheader --fai {input.fai} {input.vcf} | bcftools view -Ob > {output}
        """
