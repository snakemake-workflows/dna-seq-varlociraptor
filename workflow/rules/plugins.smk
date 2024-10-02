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


use rule tabix_known_variants as tabix_revel_scores with:
    input:
        "resources/revel_scores.tsv.gz",
    output:
        "resources/revel_scores.tsv.gz.tbi",
    params:
        get_tabix_revel_params(),
    log:
        "logs/tabix/revel.log",
