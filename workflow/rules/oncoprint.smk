rule build_oncoprint_table:
    input:
        bcf=expand("merged-calls/{group}.{{event}}.fdr-controlled.bcf", group=samples["group"].unique())
    output:
        "plots/oncoprint/{event}.tsv"
    conda:
        "../envs/oncoprinttable.yaml"
    script:
        "../scripts/build_oncoprint_matrix.py"

rule plot_oncoprint:
    input:
        "plots/oncoprint/{event}.tsv"
    output:
        report("plots/oncoprint/{event}.pdf", category="Oncoprint", caption="../report/oncoprint.rst")
    conda:
        "../envs/oncoprint.yaml"
    script:
        "../scripts/oncoprint.R"
    
