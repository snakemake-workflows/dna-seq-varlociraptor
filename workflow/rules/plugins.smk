rule download_revel:
    output:
        temp("resources/revel_scores.zip"),
    shell:
        "curl https://rothsj06.u.hpc.mssm.edu/revel-v1.3_all_chromosomes.zip -o {output}"


rule process_revel_scores:
    input:
        "resources/revel_scores.zip",
    output:
        "resources/revel_scores.tsv.gz",
    shell:
        """
        unzip -p {input} | tr "," "\t" | sed '1s/.*/#&/' | bgzip -c > results/plugins/revel_scores_tmp.tsv.gz
        zgrep -h -v ^#chr results/plugins/revel_scores_tmp.tsv.gz | awk '$3 != "." ' | sort -k1,1 -k3,3n - | cat <(zcat results/plugins/revel_scores_tmp.tsv.gz | head -n1) - | bgzip -c > {output}
        rm results/plugins/revel_scores_tmp.tsv.gz
        """
