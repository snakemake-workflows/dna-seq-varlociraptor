rule gather_germline_calls:
    input:
        calls="results/calls/{group}.freebayes.{scatteritem}.bcf",
        idx="results/calls/{group}.freebayes.{scatteritem}.bcf.csi",
    output:
        pipe("results/germline-snvs/{group}.germline_snv_candidates.bcf"),
    log:
        "logs/germline-snvs/gather-calls/{group}.log",
    params:
        extra="-a",
    wrapper:
        "v2.3.2/bio/bcftools/concat"


rule control_fdr_germline_snvs:
    input:
        "results/germline-snvs/{group}.germline_snv_candidates.bcf",
    output:
        "results/germline-snvs/{group}.bcf",
    params:
        events=germline_events,
    log:
        "logs/germline-snvs/fdr-control/{group}.log",
    conda:
        "../envs/varlociraptor.yaml"
    shell:
        "varlociraptor filter-calls control-fdr {input} --events {params.events} --mode local-smart --var snv --fdr 0.05 > {output} 2> {log}"
