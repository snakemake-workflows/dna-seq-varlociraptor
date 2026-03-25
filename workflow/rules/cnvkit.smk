# These rules are based on the existing workflow of Dawid Krzeciesa (@dawidkrzeciesa):
# https://github.com/dawidkrzeciesa/dna-seq-cnv/
# Extended by David Lähnemann (@dlaehnemann) in a feature branch:
# https://github.com/dawidkrzeciesa/dna-seq-cnv/tree/crc_exome
# Rules have been selected, edited and extended to fit the surrounding pipeline

import math
from os import path


rule cnvkit_access:
    input:
        fasta=genome,
    output:
        "results/cnvkit/access-mappable.bed",
    log:
        "logs/cnvkit/acces/access.log",
    conda:
        "../envs/cnvkit.yaml"
    params:
        access_param=lookup(dpath="params/cnvkit/access_param", within=config),
    shell:
        "(cnvkit.py access {input} {params.access_param} -o {output}) 2>{log}"


rule download_annotation_gtf:
    output:
        "resources/cnvkit/annotation.gtf",
    log:
        "logs/cnvkit/download_annotation_gff3.log",
    cache: "omit-software"  # save space and time with between workflow caching (see docs)
    params:
        species=lookup(within=config, dpath="ref/species"),
        build=lookup(within=config, dpath="ref/build"),
        release=lookup(within=config, dpath="ref/release"),
        flavor="",  # optional, e.g. chr_patch_hapl_scaff, see Ensembl FTP.
        # branch="plants",  # optional: specify branch
    wrapper:
        "v7.0.0/bio/reference/ensembl-annotation"


rule ensure_gene_name_for_all_entries:
    input:
        "resources/cnvkit/annotation.gtf",
    output:
        "resources/cnvkit/annotation.gene_name_for_all_entries.gtf",
    log:
        "logs/cnvkit/annotation.gene_name_for_all_entries.gtf",
    threads: 4
    params:
        extra='--itsv --implicit-tsv-header --pass-comments --otsv --headerless-tsv-output put \'if (!strmatch($9,"gene_name")) {$9 = sub($9, "gene_id (\\"ENSG[^\\"]+\\");", "gene_id \\1; gene_name \\1;")}\'',
    wrapper:
        "v7.1.0/utils/miller"


# do the conversion as suggested here:
# https://github.com/etal/cnvkit/issues/311#issuecomment-367500456
rule convert_gtf_to_bed:
    input:
        gtf="resources/cnvkit/annotation.gene_name_for_all_entries.gtf",
    output:
        bed="resources/cnvkit/annotation.bed",
    log:
        "logs/cnvkit/annotation.gff3_to_bed.log",
    conda:
        "../envs/cnvkit.yaml"
    resources:
        mem_mb=5000,
    shell:
        "( skg_convert.py {input.gtf} "
        "    --from gff "
        "    --to bed4 "
        "    --gff-tag gene_name "
        "    --gff-type gene "
        "    --flatten "
        "    --output {output.bed} "
        ") >{log} 2>&1"


# TODO: Produces two unused files that would require to introduce another wildcard,
# which bloats the file paths and outsources the normal-sample resolution to the
# calling rule. Maybe split rule into individual cnvkit calls.
rule cnvkit_batch:
    input:
        tumor_bam=get_cnvkit_batch_input,
        normal_bam=lambda w: get_cnvkit_batch_input(w, sample_type="normal"),
        tumor_bai=lambda w: get_cnvkit_batch_input(w, ext="bai"),
        normal_bai=lambda w: get_cnvkit_batch_input(w, sample_type="normal", ext="bai"),
        fasta=genome,
        targets=lookup(dpath="params/cnvkit/target_bed", within=config),
        access=rules.cnvkit_access.output,
        bed="resources/cnvkit/annotation.bed",
    output:
        callns="results/cnvkit/batch/{group}/{sample}/{sample}.call.cns",
        bintest="results/cnvkit/batch/{group}/{sample}/{sample}.bintest.cns",
        cnr="results/cnvkit/batch/{group}/{sample}/{sample}.cnr",
        anticnn="results/cnvkit/batch/{group}/{sample}/{sample}.antitargetcoverage.cnn",
        targetcnn="results/cnvkit/batch/{group}/{sample}/{sample}.targetcoverage.cnn",
        #"results/cnvkit/batch/{group}/{sample}/{normal}.antitargetcoverage.cnn",
        #"results/cnvkit/batch/{group}/{sample}/{normal}.targetcoverage.cnn",
        diagram="results/cnvkit/batch/{group}/{sample}/{sample}-diagram.pdf",
        scatter="results/cnvkit/batch/{group}/{sample}/{sample}-scatter.png",
        antibed="results/cnvkit/batch/{group}/{sample}/coding_regions.antitarget.bed",
        targetbed="results/cnvkit/batch/{group}/{sample}/coding_regions.target.bed",
        cns="results/cnvkit/batch/{group}/{sample}/{sample}.cns",
        cnn="results/cnvkit/batch/{group}/{sample}/{sample}.cnn",
    log:
        "logs/cnvkit-batch/{group}/{sample}.log",
    group:
        "cnvkit_{group}"
    conda:
        "../envs/cnvkit.yaml"
    threads: 64
    params:
        normal=lambda wc, input: (
            f"--normal {input.normal_bam}"
            if input.normal_bam != input.tumor_bam
            else ""
        ),
        batch_param=lookup(dpath="params/cnvkit/batch_param", within=config),
        chr_sex=get_cnvkit_sex,
        targets=lambda wc, input: (
            f"--targets {input.targets}" if input.targets else "--method wgs"
        ),
        out_dir=lambda wc, output: path.dirname(output.cns),
    shell:
        "cnvkit.py batch {input.tumor_bam}"
        "  {params.normal}"
        "  {params.targets}"
        "  --fasta {input.fasta}"
        "  --output-reference {output.cnn}"
        "  --access {input.access}"
        "  --annotate {input.bed}"
        "  --output-dir {params.out_dir}"
        "  --diagram"
        "  --scatter"
        "  -p {threads}"
        "  {params.chr_sex}"
        "  {params.batch_param};"
        "  2>{log}"


rule genotype_snvs_to_vcf:
    input:
        bcf=expand(
            "results/calls/fdr-controlled/{{group}}/{event}/{{group}}.SNV.variants.bcf",
            event=lookup(dpath="params/cnvkit/joint_event", within=config),
        ),
    output:
        vcf=temp("results/cnvkit/annotate/{group}.genotyped.vcf"),
    log:
        "logs/genotype-snvs/{group}.log",
    conda:
        "../envs/varlociraptor.yaml"
    shell:
        "varlociraptor genotype < {input.bcf} | bcftools view -O v -o {output.vcf}"


rule create_allele_depth_annotation:
    input:
        vcf="results/cnvkit/annotate/{group}.genotyped.vcf",
    output:
        tsv=temp("results/cnvkit/annotate/{group}.allele_depths.tsv"),
    log:
        "logs/vembrane-ad/{group}.log",
    conda:
        "../envs/vembrane.yaml"
    shell:
        "vembrane table "
        "  --wide "
        "  --header none "
        "  'CHROM,POS,REF,ALT,"
        'for_each_sample(lambda s: ",".join([ str(int(FORMAT["DP"][s] - round(FORMAT["DP"][s] * FORMAT["AF"][s]))),str(int(round(FORMAT["DP"][s] * FORMAT["AF"][s]))) ]))'
        "' "
        "{input.vcf} "
        ">{output.tsv} "
        "2>{log} "


rule bgzip_and_tabix_allele_depths:
    input:
        tsv="results/cnvkit/annotate/{group}.allele_depths.tsv",
    output:
        gz=temp("results/cnvkit/annotate/{group}.allele_depths.tsv.gz"),
        tbi=temp("results/cnvkit/annotate/{group}.allele_depths.tsv.gz.tbi"),
    log:
        "logs/tabix-ad/{group}.log",
    conda:
        "../envs/bcftools.yaml"
    shell:
        "( bgzip {input.tsv}; " "  tabix {output.gz} -b 2 -e 2;" ") 2> {log}"


rule annotate_tumor_allele_depth:
    input:
        vcf="results/cnvkit/annotate/{group}.genotyped.vcf",
        tsv="results/cnvkit/annotate/{group}.allele_depths.tsv.gz",
        tbi="results/cnvkit/annotate/{group}.allele_depths.tsv.gz.tbi",
    output:
        vcf="results/cnvkit/annotate/{group}.genotyped.allele_depths.vcf",
    log:
        "logs/annotate-ad/{group}.log",
    conda:
        "../envs/bcftools.yaml"
    shell:
        "( bcftools annotate "
        "    --annotations {input.tsv} "
        "    --columns CHROM,POS,REF,ALT,FORMAT/AD "
        "     --header-line '##FORMAT=<ID=AD,Number=R,Type=Integer,Description=\"Allelic Depth\">' "
        "    {input.vcf} "
        "    >{output.vcf} "
        ") 2>{log}"


rule cnvkit_call:
    input:
        cns="results/cnvkit/batch/{group}/{sample}/{sample}.cns",
        vcf="results/cnvkit/annotate/{group}.genotyped.allele_depths.vcf",
    output:
        cns="results/cnvkit/call/{group}/{sample}.final.cns",
    log:
        "logs/cnvkit-call/{group}/{sample}.cns.log",
    conda:
        "../envs/cnvkit.yaml"
    threads: 1
    params:
        call_param=lookup(dpath="params/cnvkit/call_param", within=config),
        tumor_purity=get_cnvkit_purity_setting,
        chr_sex=get_cnvkit_sex,
        sample_sex=get_sample_sex,
    shell:
        "(cnvkit.py call "
        "  -v {input.vcf} "
        "  {params.chr_sex} "
        "  -x {params.sample_sex} "
        "  --drop-low-coverage "
        "  -m clonal "
        "  {params.tumor_purity} "
        "  {input.cns} "
        "  -o {output.cns} "
        "  {params.call_param} "
        ") 2>{log}"
