rule orthanq_candidate_variants:
    input:
        ...
    output:
        "results/candidate-calls/orthanq.bcf",
    log:
        "logs/orthanq-candidates.log",
    shell:
        "orthanq candidates ... --out-bcf {output} 2> {log}"


rule gather_annotated_calls:
    input:
        calls=gather.calling("results/calls/{{group}}.orthanq.{scatteritem}.bcf"),
        idx=gather.calling("results/calls/{{group}}.orthanq.{scatteritem}.bcf.csi"),
    output:
        "results/calls/{group}.hla-variants.bcf",
    log:
        "logs/gather-hla-variants/{group}.log",
    params:
        extra="-a",
    group:
        "annotation"
    wrapper:
        "v2.3.2/bio/bcftools/concat"


rule orthanq_call:
    input:
        "results/calls/{group}.hla-variants.bcf",
    output:
        "results/hla-typing/{group}.{sample}.tsv",
    log:
        "logs/orthanq/{group}.{sample}.log",
    shell:
        "orthanq call {input} --sample {wildcards.sample} --out-table {output} 2> {log}"
