rule get_genome:
    output:
        genome,
    log:
        "logs/get-genome.log",
    cache: "omit-software"
    params:
        species=config["ref"]["species"],
        datatype="dna",
        build=config["ref"]["build"],
        release=config["ref"]["release"],
        chromosome=config["ref"].get("chromosome"),
    wrapper:
        "v7.3.0/bio/reference/ensembl-sequence"


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
    cache: "omit-software"
    conda:
        "../envs/samtools.yaml"
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
    cache: "omit-software"
    params:
        species=config["ref"]["species"],
        release=config["ref"]["release"],
        build=config["ref"]["build"],
        type="all",
        chromosome=config["ref"].get("chromosome"),
    wrapper:
        "v7.5.0/bio/reference/ensembl-variation"


rule get_annotation:
    output:
        "resources/annotation.gtf",
    log:
        "logs/get_annotation.log",
    cache: "omit-software"  # save space and time with between workflow caching (see docs)
    params:
        species=config["ref"]["species"],
        release=config["ref"]["release"],
        build=config["ref"]["build"],
        flavor="",  # optional, e.g. chr_patch_hapl_scaff, see Ensembl FTP.
    wrapper:
        "v7.5.0/bio/reference/ensembl-annotation"


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
    cache: "omit-software"
    conda:
        "../envs/rbt.yaml"
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
    log:
        "logs/vep/cache.log",
    cache: "omit-software"
    params:
        species=config["ref"]["species"],
        build=config["ref"]["build"],
        release=config["ref"]["release"],
    wrapper:
        "v8.0.0/bio/vep/cache"


rule get_vep_plugins:
    output:
        directory("resources/vep/plugins"),
    log:
        "logs/vep/plugins.log",
    params:
        release=config["ref"]["release"],
    wrapper:
        "v8.0.0/bio/vep/plugins"


rule get_pangenome:
    output:
        f"{pangenome_prefix}.{{ext}}",
    log:
        "logs/pangenome/{ext}.log",
    wildcard_constraints:
        ext="hapl|gbz",
    cache: "omit-software"
    params:
        url=lambda wc: get_pangenome_url(wc.ext),
    shell:
        "curl -o {output} {params.url} 2> {log}"
