rule get_genome:
    output:
        genome,
    log:
        "logs/get-genome.log",
    params:
        species=config["ref"]["species"],
        datatype="dna",
        build=config["ref"]["build"],
        release=config["ref"]["release"],
        chromosome=config["ref"].get("chromosome"),
    cache: "omit-software"
    wrapper:
        "v2.3.2/bio/reference/ensembl-sequence"


rule genome_faidx:
    input:
        genome,
    output:
        genome_fai,
    log:
        "logs/genome-faidx.log",
    cache: "omit-software"
    wrapper:
        "v2.3.2/bio/samtools/faidx"


rule genome_dict:
    input:
        genome,
    output:
        genome_dict,
    log:
        "logs/samtools/create_dict.log",
    conda:
        "../envs/samtools.yaml"
    cache: "omit-software"
    shell:
        "samtools dict {input} > {output} 2> {log} "


rule get_known_variants:
    input:
        # use fai to annotate contig lengths for GATK BQSR
        fai=genome_fai,
    output:
        vcf="resources/variation.vcf.gz",
    log:
        "logs/get-known-variants.log",
    params:
        species=config["ref"]["species"],
        release=config["ref"]["release"],
        build=config["ref"]["build"],
        type="all",
        chromosome=config["ref"].get("chromosome"),
    cache: "omit-software"
    wrapper:
        "v3.7.0/bio/reference/ensembl-variation"


rule get_annotation:
    output:
        "resources/annotation.gtf",
    params:
        species=config["ref"]["species"],
        release=config["ref"]["release"],
        build=config["ref"]["build"],
        flavor="",  # optional, e.g. chr_patch_hapl_scaff, see Ensembl FTP.
    log:
        "logs/get_annotation.log",
    cache: "omit-software"  # save space and time with between workflow caching (see docs)
    wrapper:
        "v2.3.2/bio/reference/ensembl-annotation"


rule determine_coding_regions:
    input:
        "resources/annotation.gtf",
    output:
        "resources/coding_regions.bed.gz",
    log:
        "logs/determine_coding_regions.log",
    cache: "omit-software"  # save space and time with between workflow caching (see docs)
    conda:
        "../envs/awk.yaml"
    shell:
        # filter for `exon` entries, but unclear how to exclude pseudogene exons...
        """
        ( cat {input} | \\
          awk 'BEGIN {{ IFS = "\\t"}} {{ if ($3 == "exon") {{ print $0 }} }}' | \\
          grep 'transcript_biotype "protein_coding"' | \\
          grep 'gene_biotype "protein_coding"' | \\
          awk 'BEGIN {{ IFS = "\\t"; OFS = "\\t"}}  {{ print $1,$4-1,$5 }}' | \\
          gzip > {output} \\
        ) 2> {log}
        """


rule remove_iupac_codes:
    input:
        "resources/variation.vcf.gz",
    output:
        "resources/variation.noiupac.vcf.gz",
    log:
        "logs/fix-iupac-alleles.log",
    conda:
        "../envs/rbt.yaml"
    cache: "omit-software"
    shell:
        "(rbt vcf-fix-iupac-alleles < {input} | bcftools view -Oz > {output}) 2> {log}"


rule bwa_index:
    input:
        genome,
    output:
        idx=multiext(genome, ".amb", ".ann", ".bwt", ".pac", ".sa"),
    log:
        "logs/bwa_index.log",
    cache: True
    wrapper:
        "v2.3.2/bio/bwa/index"


rule get_vep_cache:
    output:
        directory("resources/vep/cache"),
    params:
        species=config["ref"]["species"],
        build=config["ref"]["build"],
        release=config["ref"]["release"],
    log:
        "logs/vep/cache.log",
    cache: "omit-software"
    wrapper:
        "v3.3.5/bio/vep/cache"


rule get_vep_plugins:
    output:
        directory("resources/vep/plugins"),
    params:
        release=config["ref"]["release"],
    log:
        "logs/vep/plugins.log",
    wrapper:
        "v3.3.5/bio/vep/plugins"


rule get_pangenome:
    output:
        f"{pangenome_prefix}.{{ext}}",
    params:
        url=get_pangenome_url,
    log:
        "logs/pangenome/get_reference_{ext}.log",
    cache: "omit-software"
    shell:
        "curl -o {output} {params.url} 2> {log}"


## These rules create dist- and min-indexes when index graph is provided as bgz-file
## gbz graph is available for hprc-v1.1 but mapped records file during post processing (probably because of CHM13 reference)
# rule create_pangenome_dist_index:
#     input:
#         pangenome
#     output:
#         f"{pangenome_prefix}.dist"
#     conda:
#         "../envs/vg.yaml"
#     threads: max(workflow.cores, 1)  # use all available cores
#     shell:
#         "vg index -t {threads} -j {output} {input}"
# rule create_pangenome_minimizer_index:
#     input:
#         graph=pangenome,
#         dist=f"{pangenome_prefix}.dist"
#     output:
#         f"{pangenome_prefix}.min"
#     conda:
#         "../envs/vg.yaml"
#     threads: 16
#     log:
#         "logs/pangenome/minimizer_index.log"
#     shell:
#         "vg minimizer -t {threads} -d {input.dist} -o {output} {input.graph} 2> {log}"
