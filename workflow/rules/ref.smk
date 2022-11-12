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
    cache: True
    wrapper:
        "v1.2.0/bio/reference/ensembl-sequence"


rule genome_faidx:
    input:
        genome,
    output:
        genome_fai,
    log:
        "logs/genome-faidx.log",
    cache: True
    wrapper:
        "v1.10.0/bio/samtools/faidx"


rule genome_dict:
    input:
        genome,
    output:
        genome_dict,
    log:
        "logs/samtools/create_dict.log",
    conda:
        "../envs/samtools.yaml"
    cache: True
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
    cache: True
    wrapper:
        "v1.12.0/bio/reference/ensembl-variation"


rule get_hprc_pangenome:
    output:
        "resources/hprc-pangenome.vcf.gz",
    log:
        "logs/get-hprc-pangenome.log",
    conda:
        "../envs/curl.yaml"
    cache: "omit-software"
    shell:
        "curl -L "
        "https://s3-us-west-2.amazonaws.com/human-pangenomics/pangenomes/freeze/freeze1/minigraph-cactus/hprc-v1.0-mc.grch38.vcfbub.a100k.wave.vcf.gz "
        "> {output} 2> {log}"


rule get_annotation_gz:
    output:
        "resources/annotation.gtf.gz",
    params:
        species=config["ref"]["species"],
        release=config["ref"]["release"],
        build=config["ref"]["build"],
        flavor="",  # optional, e.g. chr_patch_hapl_scaff, see Ensembl FTP.
    log:
        "logs/get_annotation.log",
    cache: True  # save space and time with between workflow caching (see docs)
    wrapper:
        "v1.5.0/bio/reference/ensembl-annotation"


rule determine_coding_regions:
    input:
        "resources/annotation.gtf.gz",
    output:
        "resources/coding_regions.bed.gz",
    log:
        "logs/determine_coding_regions.log",
    cache: True  # save space and time with between workflow caching (see docs)
    conda:
        "../envs/awk.yaml"
    shell:
        # filter for `exon` entries, but unclear how to exclude pseudogene exons...
        """
        ( zcat {input} | \\
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
    cache: True
    shell:
        "(rbt vcf-fix-iupac-alleles < {input} | bcftools view -Oz > {output}) 2> {log}"


rule vg_index:
    input:
        genome=genome,
        # Use the HPRC human pangenome if possible. Otherwise known variation from Ensembl.
        # The latter has the disadvantage that it does not contain haplotypes.
        # See https://github.com/vgteam/vg/issues/3776.
        variants="resources/hprc-pangenome.vcf.gz"
        if config["ref"]["species"] == "homo_sapiens"
        and config["ref"]["build"] == "GRCh38"
        else "resources/variation.vcf.gz",
    output:
        f"{genome}.gbz",
    log:
        "logs/vg/index.log",
    params:
        prefix=lambda w, input: os.path.splitext(input.genome)[0],
    conda:
        "../envs/vg.yaml"
    cache: True
    threads: workflow.cores
    shell:
        "vg autoindex -t {threads} --workflow giraffe "
        "-r {input.genome} -v {input.variants} -p x "
        "--prefix {params.prefix} 2> {log}"


rule get_vep_cache:
    output:
        directory("resources/vep/cache"),
    params:
        species=config["ref"]["species"],
        build=config["ref"]["build"],
        release=config["ref"]["release"],
    log:
        "logs/vep/cache.log",
    wrapper:
        "v1.12.0/bio/vep/cache"


rule get_vep_plugins:
    output:
        directory("resources/vep/plugins"),
    params:
        release=config["ref"]["release"],
    log:
        "logs/vep/plugins.log",
    wrapper:
        "v1.12.0/bio/vep/plugins"
