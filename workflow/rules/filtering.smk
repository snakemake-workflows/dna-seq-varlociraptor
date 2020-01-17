rule filter_by_annotation:
    input:
        get_annotated_bcf
    output:
        "calls/{group}.{filter}.filtered.bcf"
    params:
        filter=lambda w: config["calling"]["filter"][w.filter]
    conda:
        "../envs/snpsift.yaml"
    shell:
        "bcftools view {input} | SnpSift filter \"{params.filter}\" | bcftools view -Ob > {output}"


rule control_fdr:
    input:
        "calls/{group}.{filter}.filtered.bcf"
    output:
        "calls/{group}.{vartype}.{event}.{filter}.fdr-controlled.bcf"
    params:
        threshold=config["calling"]["fdr-control"]["threshold"],
        events=lambda wc: config["calling"]["fdr-control"]["events"][wc.event]["varlociraptor"]
    conda:
        "../envs/varlociraptor.yaml"
    shell:
        "varlociraptor filter-calls control-fdr {input} --var {wildcards.vartype} "
        "--events {params.events} --fdr {params.threshold} > {output}"


def get_merge_input(ext=".bcf"):
    def inner(wildcards):
        return expand("calls/{{group}}.{vartype}.{{event}}.{filter}.fdr-controlled{ext}",
                      ext=ext,
                      vartype=["SNV", "INS", "DEL", "MNV"],
                      filter=config["calling"]["fdr-control"]["events"][wildcards.event]["filter"])
    return inner


rule merge_calls:
    input:
        calls=get_merge_input(".bcf"),
        idx=get_merge_input(".bcf.csi")
    output:
        "merged-calls/{group}.{event}.fdr-controlled.bcf"
    params:
        "-a -Ob"
    wrapper:
        "0.37.1/bio/bcftools/concat"
