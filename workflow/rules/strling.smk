rule strling_call:
    input:
        bam="results/recal/{sample}.bam",
        bin="results/strling/extract/{sample}.bin",
        merged="results/strling/merge/merged-bounds.txt",
        ref="results/refs/genome.fasta",
    output:
        bounds="results/strling_call/{sample}-bounds.txt",
        genotype="results/strling_call/{sample}-genotype.txt",
        unplaced="results/strling_call/{sample}-unplaced.txt",
    params:
        prefix="results/strling_call/{sample}"
    log:
        "logs/strling/{sample}.call.log"
    shell:
        "workflow/scripts/strling call -f {input.ref} {input.bam} {input.bin} -o {params.prefix} -b {input.merged} 2> {log}"


rule strling_merge:
    input:
        bam=expand("results/strling/extract/{sample}.bin", sample=samples["sample_name"].values),
        ref="results/refs/genome.fasta",
    output:
        merged="results/strling/merge/merged-bounds.txt"
    params:
        prefix="results/strling/merge/merged"
    log:
        "logs/strling/merge.log"
    shell:
        "workflow/scripts/strling merge -f {input.ref} -o {params.prefix} {input.bam} 2> {log}"


rule strling_extract:
    input:
        bam="results/dedup/{sample}.bam",
        ref="results/refs/genome.fasta",
    output:
        bin="results/strling/extract/{sample}.bin",
    log:
        "logs/strling/{sample}.extract.log"
    shell:
        "workflow/scripts/strling extract -f {input.ref} {input.bam} {output.bin} 2> {log}"