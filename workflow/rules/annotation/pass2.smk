annotations = [(e, f) for e, f in config["annotations"]["vcfs"].items() if e != "activate"]


def get_annotation_pipes():
     if annotations:
         return "| " + " | ".join(
             ["SnpSift annotate -name {prefix}_ {path} /dev/stdin".format(prefix=prefix, path=path)
              for prefix, path in annotations]
         )
     else:
         return ""


def get_annotation_vcfs(idx=False):
    fmt = lambda f: f if not idx else f + ".tbi"
    return [fmt(f) for e, f in annotations]


#What about multiple ID Fields?
rule annotate_vcfs:
    input:
        bcf="calls/{prefix}.bcf",
        annotations=get_annotation_vcfs(),
        idx=get_annotation_vcfs(idx=True)
    output:
        "calls/{prefix}.db-annotated.bcf"
    params:
        extra="-Xmx4g",
        pipes=get_annotation_pipes()
    conda:
        "../../envs/snpsift.yaml"
    shell:
        "bcftools view {input.bcf} {params.pipes} | bcftools view -Ob > {output}"


rule annotate_dgidb:
    input:
        "calls/{prefix}.bcf"
    output:
        "calls/{prefix}.dgidb.bcf"
    params:
    conda:
        "../../envs/annotate_dgidb.yaml"
    resources:
        dgidb_requests=1
    shell:
        "rbt vcf-annotate-dgidb {input} > {output}"


if is_activated("annotations/dbnsfp"):
    #TODO create wrapper of this entire block
    rule dbnsfp_download:
        output:
            "resources/dbnsfp.zip"
        params:
            zip = config["annotations"]["dbnsfp"]["url"]
        shell:
            "wget {params.zip} -O {output}"


    rule dbnsfp_bgzip:
        input:
            "resources/dbnsfp.zip"
        output:
            "resources/dbnsfp.txt.gz"
        conda:
            "../../envs/htslib.yaml"
        threads: 4
        shell:
            """
            unzip {input} -d resources/dbnsfp/ &&
            (zcat resources/dbnsfp/*_variant.chr1.gz |
            head -n 1 ; zcat resources/dbnsfp/*_variant.chr* |
            grep -v '^#' ) | bgzip -@ {threads} > {output} &&
            rm -r resources/dbnsfp
            """
    
    rule annotate_dbnsfp:
        input:
            bcf="calls/{prefix}.bcf",
            db="resources/dbnsfp.txt.gz",
            idx="resources/dbnsfp.txt.gz.tbi"
        output:
            "calls/{prefix}.dbnsfp.bcf"
        params:
            extra="-Xmx4g",
            fields = ",".join(config["annotations"]["dbnsfp"]["fields"])
        conda:
            "../../envs/snpsift.yaml"
        shell:
            "bcftools view {input.bcf} | SnpSift dbnsfp -db {input.db} -f {params.fields} {params.extra} /dev/stdin | "
            "sed 's/\\(^##INFO=<ID=dbNSFP_\\w*,Number=\\)A/\\1./g' | bcftools view -Ob > {output}"


    rule create_dbnsfp_tabix_index:
        input:
            "resources/{file}.gz",
        output:
            "resources/{file}.gz.tbi"
        conda:
            "../../envs/htslib.yaml"
        shell:
            "tabix -s 1 -b 2 -e 2 {input}"
