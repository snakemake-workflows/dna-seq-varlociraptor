rule group_bcf_to_vcf:
    input:
        "results/final-calls/{group}.{event}.{calling_type}.fdr-controlled.bcf",
    output:
        temp("results/maf/{group}.{event}.{calling_type}.fdr-controlled.vcf"),
    log:
        "logs/maftools/{group}.{event}.{calling_type}.fdr-controlled.log",
    wrapper:
        "v3.8.0/bio/bcftools/view"


rule group_vcf_to_maf:
    input:
        vcf="results/maf/{group}.{event}.{calling_type}.fdr-controlled.vcf",
        ref=genome,
    output:
        maf="results/maf/{group}.{event}.{calling_type}.fdr-controlled.maf",
    log:
        "logs/maftools/{group}.{event}.{calling_type}.fdr-controlled.log",
    conda:
        "../envs/vcf2maf.yaml"
    params:
        genome_build=config["ref"]["build"],
        species=config["ref"]["species"],
    shell:
        "vcf2maf.pl --inhibit-vep "
        " --retain-info PROB_ARTIFACT,PROB_ABSENT,PROB_SOMATIC,PROB_GERMLINE_HOM,PROB_GERMLINE_HET "
        " --retain-fmt DP,AF "
        " --ncbi-build {params.genome_build} "
        " --species {params.species} "
        " --vcf-tumor-id tumor --vcf-normal-id panel_of_normals "
        " --input-vcf {input.vcf} --output-maf {output.maf} --ref-fasta {input.ref} "
        " --tumor-id {wildcards.group}_tumor --normal-id {wildcards.group}_panel_of_normal "
        " 2>{log}"
