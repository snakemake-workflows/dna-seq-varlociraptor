rule filter_candidates_by_annotation:
    input:
        "results/candidate-calls/{group}.{caller}.{scatteritem}.annotated.bcf"
    output:
        "results/candidate-calls/{group}.{caller}.{scatteritem}.filtered.bcf"
    log:
        "logs/filter-calls/annotation/{group}.{caller}.{scatteritem}.log"
    params:
        filter=lambda w: config["calling"]["filter"]["candidates"]
    conda:
        "../envs/vembrane.yaml"
    shell:
        "vembrane filter {params.filter:q} {input} --output-fmt bcf --output {output} &> {log}"


rule filter_by_annotation:
    input:
        get_annotated_bcf
    output:
        "results/calls/{group}.{filter}.{scatteritem}.filtered_ann.bcf"
    log:
        "logs/filter-calls/annotation/{group}.{filter}.{scatteritem}.log"
    params:
        filter=lambda w: config["calling"]["filter"][w.filter]
    conda:
        "../envs/vembrane.yaml"
    shell:
        "vembrane filter {params.filter:q} {input} --output-fmt bcf --output {output} &> {log}"


rule filter_odds:
    input:
        "results/calls/{group}.{filter}.{scatteritem}.filtered_ann.bcf"
    output:
        "results/calls/{group}.{event}.{filter}.{scatteritem}.filtered_odds.bcf"
    params:
        events=lambda wc: config["calling"]["fdr-control"]["events"][wc.event]["varlociraptor"]
    log:
        "logs/filter-calls/posterior_odds/{group}.{event}.{filter}.{scatteritem}.log"
    conda:
        "../envs/varlociraptor.yaml"
    shell:
        "varlociraptor filter-calls posterior-odds --events {params.events} --odds barely < {input} > {output} 2> {log}"


rule gather_calls:
    input:
        calls=gather.calling("results/calls/{{group}}.{{event}}.{{filter}}.{scatteritem}.filtered_odds.bcf"),
        idx=gather.calling("results/calls/{{group}}.{{event}}.{{filter}}.{scatteritem}.filtered_odds.bcf.csi"),
    output:
        "results/calls/{group}.{event}.{filter}.filtered_odds.bcf"
    log:
        "logs/gather-calls/{group}.{event}.{filter}.log"
    params:
        "-a -Ob"
    wrapper:
        "0.67.0/bio/bcftools/concat"


rule control_fdr:
    input:
        ("results/calls/{group}.{event}.{filter}.filtered_odds.bcf"
         if not is_activated("benchmarking") else "results/calls/{group}.bcf")
    output:
        "results/calls/{group}.{vartype}.{event}.{filter}.fdr-controlled.bcf"
    log:
        "logs/control-fdr/{group}.{vartype}.{event}.{filter}.log"
    params:
        query=get_fdr_control_params
    conda:
        "../envs/varlociraptor.yaml"
    shell:
        "varlociraptor filter-calls control-fdr {input} --var {wildcards.vartype} "
        "--events {params.query[events]} --fdr {params.query[threshold]} > {output} 2> {log}"


rule merge_calls:
    input:
        calls=get_merge_calls_input(".bcf"),
        idx=get_merge_calls_input(".bcf.csi")
    output:
        "results/merged-calls/{group}.{event}.fdr-controlled.bcf"
    log:
        "logs/merge-calls/{group}.{event}.log"
    params:
        "-a -Ob"
    wrapper:
        "0.59.2/bio/bcftools/concat"
