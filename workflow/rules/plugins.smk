rule download_revel:
    output:
        temp("resources/revel_scores.zip"),
    log:
        "logs/vep_plugins/download_revel.log",
    conda:
        "../envs/curl.yaml"
    shell:
        "curl https://zenodo.org/records/7072866/files/revel-v1.3_all_chromosomes.zip -o {output} &> {log}"


rule process_revel_scores:
    input:
        "resources/revel_scores.zip",
    output:
        "resources/revel_scores.tsv.gz",
    params:
        build=config["ref"]["build"],
    log:
        "logs/vep_plugins/process_revel_scores.log",
    conda:
        "../envs/htslib.yaml"
    shell:
        """
        tmpfile=$(mktemp {resources.tmpdir}/revel_scores.XXXXXX)
        unzip -p {input} | tr "," "\t" | sed '1s/.*/#&/' | bgzip -c > $tmpfile
        if [ "{params.build}" == "GRCh38" ] ; then
            zgrep -h -v ^#chr $tmpfile | awk '$3 != "." ' | sort -k1,1 -k3,3n - | cat <(zcat $tmpfile | head -n1) - | bgzip -c > {output}
        elif [ "{params.build}" == "GRCh37" ] ; then
            cat $tmpfile > {output}
        else
            echo "Annotation of REVEL scores only supported for GRCh37 or GRCh38" > {log}
            exit 125
        fi
        """


rule download_cadd_scores_for_vep:
    output:
        cadd="resources/cadd/{build}/{cadd_version}/{variant_type}.tsv.gz",
    log:
        "logs/cadd/{build}/{cadd_version}/{variant_type}.log",
    params:
        file_name=lambda wc: (
            "whole_genome_SNVs"
            if wc.variant_type == "snv"
            else (
                "gnomad.genomes.r4.0.indel"
                if wc.variant_type == "indels"
                else "unknown_variant_type_choose_snvs_or_indels"
            )
        ),
    cache: "omit-software"
    conda:
        "../envs/download_cadd.yaml"
    shell:
        "( wget --retry-connrefused --waitretry=10 --tries=10 --continue "
        '    "https://kircherlab.bihealth.org/download/CADD/{wildcards.cadd_version}/{wildcards.build}/{params.file_name}.tsv.gz" '
        "    -O - | "
        "   gzip -d | "
        "   cut -f1-4,6 | "
        "   bgzip -l1 -c "
        "   > {output.cadd} "
        ") 2>{log} "
