rule download_snpeff_db:
    output:
        directory("resources/snpeff/{ref}")
    log:
        "logs/download-snpeff-db/{ref}.log"
    params:
        db_dir=lambda _, output: str(Path(output[0]).parent.resolve()),
        ref="{ref}"
    conda:
        "../envs/snpeff.yaml"
    cache: True
    shell:
        "snpEff download -dataDir {params.db_dir} {params.ref} 2> {log}"

rule snpeff:
    input:
        calls="results/calls/{group}.bcf",
        db="resources/snpeff/{build}.{snpeff_release}".format(**config["ref"])
    output:
        calls="results/calls/{group}.annotated.bcf",
        stats="results/snpeff/{group}.html",
        csvstats="results/snpeff/{group}.csv"
    log:
        "logs/snpeff/{group}.log"
    params:
        reference="{build}.{snpeff_release}".format(**config["ref"]),
        data_dir=lambda _, input: Path(input.db).parent.resolve(),
        extra="-Xmx4g -nodownload"
    resources:
        mem_mb=4000
    wrapper:
        "0.50.4/bio/snpeff"

# TODO What about multiple ID Fields?
rule annotate_vcfs:
    threads:
        100
    input:
        bcf="results/calls/{prefix}.bcf",
        csi="results/calls/{prefix}.bcf.csi",
        annotations=get_annotation_vcfs(),
        idx=get_annotation_vcfs(idx=True)
    output:
        "results/calls/{prefix}.db-annotated.bcf"
    log:
        "logs/annotate-vcfs/{prefix}.log"
    params:
        extra="-Xmx4g",
        pipes=get_annotation_pipes
    conda:
        "../envs/snpsift.yaml"
    shell:
        "(python ../workflow/scripts/parallel_vcf.py --threads {threads} {input.bcf} '{params.pipes}' | bcftools view --threads {threads} -Ob > {output}) 2> {log}"


rule annotate_dgidb:
    threads:
        4
    input:
        "results/calls/{prefix}.bcf"
    output:
        "results/calls/{prefix}.dgidb.bcf"
    log:
        "logs/annotate-dgidb/{prefix}.log"
    conda:
        "../envs/annotate_dgidb.yaml"
    resources:
        dgidb_requests=1
    shell:
        "rbt vcf-annotate-dgidb {input} > {output} 2> {log}"


if is_activated("annotations/dbnsfp"):
    #TODO create wrapper of this entire block
    rule download_dbnsfp:
        output:
            "resources/dbnsfp.zip"
        log:
            "logs/download-dbnsfp.log"
        params:
            zip = config["annotations"]["dbnsfp"]["url"]
        shell:
            "wget {params.zip} -O {output} 2> {log}"


    rule extract_dbnsfp:
        threads:
             4
        input:
            "resources/dbnsfp.zip"
        output:
            "resources/dbnsfp.txt.gz"
        log:
            "logs/extract-dbnsfp.log"
        conda:
            "../envs/htslib.yaml"
        cache: True
        shell:
            "cat <(unzip -p {input} '*_variant.chr1.gz' | zcat | "
            "head -n 1) <(unzip -p {input} '*_variant.chr*' | "
            "zcat | grep -v '^#') | bgzip -l9 -@ {threads} > {output} 2> {log}"
    
    rule annotate_dbnsfp:
        threads:
            4
        input:
            bcf="results/calls/{prefix}.bcf",
            db="resources/dbnsfp.txt.gz",
            idx="resources/dbnsfp.txt.gz.tbi"
        output:
            "results/calls/{prefix}.dbnsfp.bcf"
        log:
            "logs/annotate-dbnsfp/{prefix}.log"
        params:
            extra="-Xmx4g",
            fields = ",".join(config["annotations"]["dbnsfp"]["fields"])
        conda:
            "../envs/snpsift.yaml"
        shell:
            "(bcftools view --threads {threads} {input.bcf} | SnpSift dbnsfp -db {input.db} -f {params.fields} {params.extra} /dev/stdin | "
            "sed 's/\\(^##INFO=<ID=dbNSFP_\\w*,Number=\\)A/\\1./g' | bcftools view -Ob --threads {threads} > {output}) 2> {log}"
