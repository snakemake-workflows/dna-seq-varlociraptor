rule filter_candidates_by_annotation:
    input:
        "results/candidate-calls/{group}.{caller}.{scatteritem}.annotated.bcf",
    output:
        "results/candidate-calls/{group}.{caller}.{scatteritem}.filtered.bcf",
    log:
        "logs/filter-calls/annotation/{group}.{caller}.{scatteritem}.log",
    params:
        filter=lambda w: config["calling"]["filter"]["candidates"],
    conda:
        "../envs/vembrane.yaml"
    shell:
        "(bcftools norm -Ou --do-not-normalize --multiallelics -any {input} | "
        "vembrane filter {params.filter:q} --output-fmt bcf --output {output}) &> {log}"


rule filter_by_annotation:
    input:
        get_annotated_bcf,
    output:
        "results/calls/{group}.{event}.{scatteritem}.filtered_ann.bcf",
    log:
        "logs/filter-calls/annotation/{group}.{event}.{scatteritem}.log",
    params:
        filter=get_annotation_filter,
    conda:
        "../envs/vembrane.yaml"
    shell:
        "vembrane filter {params.filter:q} {input} --output-fmt bcf --output {output} &> {log}"


rule filter_odds:
    input:
        "results/calls/{group}.{event}.{scatteritem}.filtered_ann.bcf",
    output:
        "results/calls/{group}.{event}.{scatteritem}.filtered_odds.bcf",
    params:
        events=lambda wc: config["calling"]["fdr-control"]["events"][wc.event][
            "varlociraptor"
        ],
    log:
        "logs/filter-calls/posterior_odds/{group}.{event}.{scatteritem}.log",
    conda:
        "../envs/varlociraptor.yaml"
    shell:
        "varlociraptor filter-calls posterior-odds --events {params.events} --odds barely < {input} > {output} 2> {log}"


rule gather_calls:
    input:
        calls=get_gather_calls_input(),
        idx=get_gather_calls_input(ext="bcf.csi"),
    output:
        "results/calls/{group}.{event}.filtered_{by}.bcf",
    log:
        "logs/gather-calls/{group}.{event}.filtered_{by}.log",
    params:
        "-a -Ob",
    wrapper:
        "v1.2.0/bio/bcftools/concat"


rule control_fdr:
    input:
        get_control_fdr_input,
    output:
        "results/calls/{group}.{vartype}.{event}.fdr-controlled.bcf",
    log:
        "logs/control-fdr/{group}.{vartype}.{event}.log",
    params:
        query=get_fdr_control_params,
    conda:
        "../envs/varlociraptor.yaml"
    shell:
        "varlociraptor filter-calls control-fdr {input} {params.query[local]} --var {wildcards.vartype} "
        "--events {params.query[events]} --fdr {params.query[threshold]} > {output} 2> {log}"


rule merge_calls:
    input:
        calls=get_merge_calls_input("bcf"),
        idx=get_merge_calls_input("bcf.csi"),
    output:
        "results/final-calls/{group}.{event}.fdr-controlled.bcf",
    log:
        "logs/merge-calls/{group}.{event}.log",
    params:
        "-a -Ob",
    wrapper:
        "v1.2.0/bio/bcftools/concat"
