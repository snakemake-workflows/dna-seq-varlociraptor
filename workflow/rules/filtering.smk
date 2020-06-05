rule filter_by_annotation:
    input:
        get_annotated_bcf
    output:
        "results/calls/{group}.{filter}.filtered.bcf"
    log:
        "logs/filter-calls/{group}.{filter}.log"
    params:
        filter=lambda w: config["calling"]["filter"][w.filter]
    conda:
        "../envs/snpsift.yaml"
    shell:
        "(bcftools view {input} | SnpSift filter \"{params.filter}\" | bcftools view -Ob > {output}) 2> {log}"


rule control_fdr:
    input:
        "results/calls/{group}.{filter}.filtered.bcf"
    output:
        "results/calls/{group}.{vartype}.{event}.{filter}.fdr-controlled.bcf"
    log:
        "logs/control-fdr/{group}.{vartype}.{event}.{filter}.log"
    params:
        threshold=config["calling"]["fdr-control"]["threshold"],
        events=lambda wc: config["calling"]["fdr-control"]["events"][wc.event]["varlociraptor"]
    conda:
        "../envs/varlociraptor.yaml"
    shell:
        "varlociraptor filter-calls control-fdr {input} --var {wildcards.vartype} "
        "--events {params.events} --fdr {params.threshold} > {output} 2> {log}"


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
