rule group_bcf_to_vcf:
    input:
        "results/final-calls/{group}.{event}.fdr-controlled.bcf",
    output:
        temp("results/final-calls/{group}.{event}.fdr-controlled.vcf"),
    log:
        "logs/final-calls/{group}.{event}.fdr-controlled.log",
    params:
        extra="",
    wrapper:
        "v3.8.0/bio/bcftools/view"

rule group_vcf_to_maf:
    input:
        vcf="results/final-calls/{group}.{event}.fdr-controlled.vcf",
        ref=dna_seq_varlociraptor.genome,
    output:
        maf="results/final-calls/{group}.{event}.fdr-controlled.maf",
    log:
        "logs/final-calls/{group}.{event}.fdr-controlled.log"
    conda:
        "envs/vcf2maf.yaml"
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
