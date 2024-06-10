rule filter_candidates_by_annotation:
    input:
        bcf="results/candidate-calls/{group}.{caller}.{scatteritem}.annotated.bcf",
        aux=get_candidate_filter_aux_files(),
    output:
        "results/candidate-calls/{group}.{caller}.{scatteritem}.filtered.bcf",
    log:
        "logs/filter-calls/annotation/{group}.{caller}.{scatteritem}.log",
    params:
        filter=get_candidate_filter_expression,
        aux=get_candidate_filter_aux(),
    conda:
        "../envs/vembrane.yaml"
    shell:
        "(bcftools norm -Ou --do-not-normalize --multiallelics -any {input} | "
        'vembrane filter {params.aux} "{params.filter}" | bcftools sort -Ob > {output}) &> {log}'


rule filter_by_annotation:
    input:
        bcf=get_annotated_bcf,
        aux=get_annotation_filter_aux_files,
    output:
        "results/calls/{group}.{event}.{calling_type}.{scatteritem}.filtered_ann.bcf",
    log:
        "logs/filter-calls/annotation/{group}.{event}.{calling_type}.{scatteritem}.log",
    params:
        filter=get_annotation_filter_expression,
        aux=get_annotation_filter_aux,
    conda:
        "../envs/vembrane.yaml"
    shell:
        'vembrane filter {params.aux} "{params.filter}" {input.bcf} --output-fmt bcf --output {output} &> {log}'


rule filter_odds:
    input:
        "results/calls/{group}.{event}.{calling_type}.{scatteritem}.filtered_ann.bcf",
    output:
        "results/calls/{group}.{event}.{calling_type}.{scatteritem}.filtered_odds.bcf",
    params:
        events=lambda wc: config["calling"]["fdr-control"]["events"][wc.event][
            "varlociraptor"
        ],
    log:
        "logs/filter-calls/posterior_odds/{group}.{event}.{calling_type}.{scatteritem}.log",
    conda:
        "../envs/varlociraptor.yaml"
    shell:
        "varlociraptor filter-calls posterior-odds --events {params.events} --odds barely < {input} > {output} 2> {log}"


rule gather_calls:
    input:
        calls=get_gather_calls_input(),
        idx=get_gather_calls_input(ext="bcf.csi"),
    output:
        "results/calls/{group}.{event}.{calling_type}.filtered_{by}.bcf",
    log:
        "logs/gather-calls/{group}.{event}.{calling_type}.filtered_{by}.log",
    params:
        extra="-a",
    wrapper:
        "v2.3.2/bio/bcftools/concat"


rule control_fdr:
    input:
        get_control_fdr_input,
    output:
        "results/calls/{group}.{vartype}.{event}.{calling_type}.fdr-controlled.bcf",
    log:
        "logs/control-fdr/{group}.{vartype}.{event}.{calling_type}.log",
    params:
        query=get_fdr_control_params,
    conda:
        "../envs/varlociraptor.yaml"
    shell:
        "varlociraptor filter-calls control-fdr {input} {params.query[mode]} --var {wildcards.vartype} "
        "--events {params.query[events]} --fdr {params.query[threshold]} > {output} 2> {log}"


rule merge_calls:
    input:
        calls=get_merge_calls_input("bcf"),
        idx=get_merge_calls_input("bcf.csi"),
    output:
        "results/final-calls/{group}.{event}.{calling_type}.fdr-controlled.bcf",
    log:
        "logs/merge-calls/{group}.{event}.{calling_type}.log",
    params:
        extra="-a",
    wrapper:
        "v2.3.2/bio/bcftools/concat"


rule convert_phred_scores:
    input:
        "results/final-calls/{group}.{event}.{calling_type}.fdr-controlled.bcf",
    output:
        "results/final-calls/{group}.{event}.{calling_type}.fdr-controlled.normal-probs.bcf",
    log:
        "logs/convert-phred-scores/{group}.{event}.{calling_type}.log",
    conda:
        "../envs/varlociraptor.yaml"
    shell:
        "varlociraptor decode-phred < {input} > {output} 2> {log}"
