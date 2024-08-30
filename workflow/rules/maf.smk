rule group_bcf_to_vcf:
    input:
        "results/final-calls/{group}.{event}.{calling_type}.fdr-controlled.bcf",
    output:
        temp("results/maf/{group}.{event}.{calling_type}.fdr-controlled.vcf"),
    log:
        "logs/maf/{group}.{event}.{calling_type}.fdr-controlled.log",
    wrapper:
        "v3.8.0/bio/bcftools/view"


rule group_vcf_to_maf:
    input:
        vcf="results/maf/{group}.{event}.{calling_type}.fdr-controlled.vcf",
        ref=genome,
        scenario="results/scenarios/{group}.yaml",  # needed for determining event probability INFO fields
    output:
        maf="results/maf/{group}.{event}.{calling_type}.fdr-controlled.maf",
    log:
        "logs/maf/{group}.{event}.{calling_type}.fdr-controlled.log",
    conda:
        "../envs/vcf2maf.yaml"
    params:
        genome_build=lookup(dpath="ref/build", within=config),
        species=lookup(dpath="ref/species", within=config),
        ann_fields=lambda wc, input: ",".join(get_annotation_fields_for_tables(wc)),
        format_fields=lambda wc, input: ",".join(get_format_fields_for_tables(wc)),
        info_fields=lambda wc, input: ",".join(
            get_info_prob_fields_for_tables(wc, input)
            + get_info_fusion_fields_for_tables(wc)
        ),
        vcf_primary_alias=lookup(
            dpath="maf/primary_alias", within=config, default="tumor"
        ),
        vcf_control_alias_option=(
            f'--vcf-normal-id {lookup(dpath= "maf/control_alias", within= config, default= "")}'
            if lookup(dpath="maf/control_alias", within=config, default=False)
            else ""
        ),
        normal_id=lambda wc: (
            f'--normal-id {wc.group}_{lookup(dpath= "maf/control_alias", within= config, default= "")}'
            if lookup(dpath="maf/control_alias", within=config, default=False)
            else ""
        ),
    shell:
        "vcf2maf.pl --inhibit-vep "
        " --retain-ann {params.ann_fields} "
        " --retain-fmt {params.format_fields} "
        " --retain-info {params.info_fields} "
        " --ncbi-build {params.genome_build} "
        " --species {params.species} "
        " --vcf-tumor-id {params.vcf_primary_alias} "
        " {params.vcf_control_alias_option} "
        " --input-vcf {input.vcf} "
        " --output-maf {output.maf} "
        " --ref-fasta {input.ref} "
        " --tumor-id {wildcards.group}_{params.vcf_primary_alias} "
        " {params.normal_id} "
        " 2>{log};"
        '! (grep -v "WARNING: No genotype column for NORMAL in VCF!" {log} | '
        '   grep "WARNING: No genotype column for" ) '
