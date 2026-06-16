rule enrichment_background_genes:
    input:
        coverage="results/coverage/{group}.csv",
    output:
        "results/enrichment/{group}.background_genes.txt",
    params:
        min_cov=config.get("enrichment", {}).get("min_gene_coverage", 20),
    log:
        "logs/enrichment/{group}.background_genes.log",
    conda:
        "../envs/pandas.yaml"
    script:
        "../scripts/enrichment_background_genes.py"


rule enrichment_analysis:
    input:
        gene_list="results/enrichment/{group}.gene_list.txt",
        background_genes="results/enrichment/{group}.background_genes.txt",
    output:
        enrichr_dir=directory("results/enrichment/{group}/{library}"),
    params:
        gene_sets=lambda wc: [wc.library],
        extra={"cutoff": 1},
    log:
        wrapper="logs/enrichment/{group}.{library}.wrapper.log",
        gseapy="logs/enrichment/{group}.{library}.log",
    threads: 1
    wrapper:
        "090be4b6c569fbadabd55c421eff58b277783ef6/bio/gseapy/gsea"


rule enrichment_gene_list:
    input:
        expand(
            "results/tables/{{group}}/{{group}}.{event}.{calling_type}.fdr-controlled.tsv",
            event=config["calling"]["fdr-control"]["events"],
            calling_type=["variants"],
        ),
    output:
        "results/enrichment/{group}.gene_list.txt",
    log:
        "logs/enrichment/{group}.gene_list.log",
    conda:
        "../envs/pandas.yaml"
    script:
        "../scripts/enrichment_gene_list.py"
