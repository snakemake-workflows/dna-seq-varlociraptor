import os

annotations = [(e, os.path.join("results", f)) for e, f in config["annotations"]["vcfs"].items() if e != "activate"]

def get_annotation_pipes(wildcards, input):
     if annotations:
         return "| " + " | ".join(
             ["SnpSift annotate -name {prefix}_ {path} /dev/stdin".format(prefix=prefix, path=path)
              for (prefix, _), path in zip(annotations, input.annotations)]
         )
     else:
         return ""


def get_annotation_vcfs(idx=False):
    fmt = lambda f: f if not idx else f + ".tbi"
    return [fmt(f) for _, f in annotations]


#What about multiple ID Fields?
rule annotate_vcfs:
    threads:
        4
    input:
        bcf="results/calls/{prefix}.bcf",
        annotations=get_annotation_vcfs(),
        idx=get_annotation_vcfs(idx=True)
    output:
        "results/calls/{prefix}.db-annotated.bcf"
    params:
        extra="-Xmx4g",
        pipes=get_annotation_pipes
    conda:
        "../../envs/snpsift.yaml"
    shell:
        "bcftools view --threads {threads} {input.bcf} {params.pipes} | bcftools view --threads {threads} -Ob > {output}"


rule annotate_dgidb:
    threads:
        4
    input:
        "results/calls/{prefix}.bcf"
    output:
        "results/calls/{prefix}.dgidb.bcf"
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
        threads:
             4
        input:
            "resources/dbnsfp.zip"
        output:
            "resources/dbnsfp.txt.gz"
        conda:
            "../../envs/htslib.yaml"
        shell:
            """
            (unzip -p {input} "*_variant.chr1.gz" | zcat |
            head -n 1 ; unzip -p {input} "*_variant.chr*" |
            zcat | grep -v '^#' ) | bgzip -l9 -@ {threads} > {output}
            """
    
    rule annotate_dbnsfp:
        threads:
            4
        input:
            bcf="results/calls/{prefix}.bcf",
            db="resources/dbnsfp.txt.gz",
            idx="resources/dbnsfp.txt.gz.tbi"
        output:
            "results/calls/{prefix}.dbnsfp.bcf"
        params:
            extra="-Xmx4g",
            fields = ",".join(config["annotations"]["dbnsfp"]["fields"])
        conda:
            "../../envs/snpsift.yaml"
        shell:
            "bcftools view --threads {threads} {input.bcf} | SnpSift dbnsfp -db {input.db} -f {params.fields} {params.extra} /dev/stdin | "
            "sed 's/\\(^##INFO=<ID=dbNSFP_\\w*,Number=\\)A/\\1./g' | bcftools view -Ob --threads {threads} > {output}"


    rule create_dbnsfp_tabix_index:
        input:
            "resources/{file}.gz",
        output:
            "resources/{file}.gz.tbi"
        conda:
            "../../envs/htslib.yaml"
        shell:
            "tabix -s 1 -b 2 -e 2 {input}"
