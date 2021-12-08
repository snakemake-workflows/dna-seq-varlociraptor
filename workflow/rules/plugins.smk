rule download_revel:
    output:
        temp("resources/revel_scores.zip"),
    log:
        "logs/download_revel.log"
    shell:
        "curl https://rothsj06.u.hpc.mssm.edu/revel-v1.3_all_chromosomes.zip -o {output} &> {log}"


rule process_revel_scores:
    input:
        "resources/revel_scores.zip",
    output:
        "resources/{ref}_revel_scores.tsv.gz",
    shell:
        """
        tmpfile=$(mktemp {resources.tmpdir}/revel_scores.XXXXXX)
        unzip -p {input} | tr "," "\t" | sed '1s/.*/#&/' | bgzip -c > $tmpfile
        if [ "{wildcards.ref}" == "GRCh38" ] ; then
            zgrep -h -v ^#chr $tmpfile | awk '$3 != "." ' | sort -k1,1 -k3,3n - | cat <(zcat $tmpfile | head -n1) - | bgzip -c > {output}
        elif [ "{wildcards.ref}" == "GRCh37" ] ; then
            cat $tmpfile > {output}
        else
            echo "Annotation of REVEL scores only supported for GRCh37 or GRCh38"
            exit 125
        fi
        """
