import io


rule render_scenario:
    input:
        local(config["calling"]["scenario"]),
    output:
        report(
            "results/scenarios/{group}.yaml",
            caption="../report/scenario.rst",
            category="Variant calling scenarios",
        ),
    log:
        "logs/render-scenario/{group}.log",
    params:
        samples=samples,
    conda:
        None
    script:
        "../scripts/render-scenario.py"


rule varlociraptor_alignment_properties:
    input:
        ref="resources/genome.fasta",
        ref_idx="resources/genome.fasta.fai",
        bam="results/recal/{sample}.sorted.bam",
    output:
        "results/alignment-properties/{group}/{sample}.json",
    log:
        "logs/varlociraptor/estimate-alignment-properties/{group}/{sample}.log",
    conda:
        "../envs/varlociraptor.yaml"
    shell:
        "varlociraptor estimate alignment-properties {input.ref} --bam {input.bam} > {output} 2> {log}"


rule varlociraptor_preprocess:
    input:
        ref="resources/genome.fasta",
        ref_idx="resources/genome.fasta.fai",
        candidates=get_candidate_calls(),
        bam="results/recal/{sample}.sorted.bam",
        bai="results/recal/{sample}.sorted.bai",
        alignment_props="results/alignment-properties/{group}/{sample}.json",
    output:
        "results/observations/{group}/{sample}.{caller}.{scatteritem}.bcf",
    params:
        extra=config["params"]["varlociraptor"]["preprocess"],
    log:
        "logs/varlociraptor/preprocess/{group}/{sample}.{caller}.{scatteritem}.log",
    benchmark:
        "benchmarks/varlociraptor/preprocess/{group}/{sample}.{caller}.{scatteritem}.tsv"
    conda:
        "../envs/varlociraptor.yaml"
    shell:
        "varlociraptor preprocess variants --candidates {input.candidates} {params.extra} "
        "--alignment-properties {input.alignment_props} {input.ref} --bam {input.bam} --output {output} "
        "2> {log}"


rule varlociraptor_call:
    input:
        obs=get_group_observations,
        scenario="results/scenarios/{group}.yaml",
    output:
        temp("results/calls/{group}.{caller}.{scatteritem}.bcf"),
    log:
        "logs/varlociraptor/call/{group}.{caller}.{scatteritem}.log",
    params:
        obs=lambda w, input: [
            "{}={}".format(s, f) for s, f in zip(get_group_aliases(w), input.obs)
        ],
        extra=config["params"]["varlociraptor"]["call"],
        postprocess=">"
        if not config["calling"].get("infer_genotypes")
        else "| varlociraptor genotype >",
    conda:
        "../envs/varlociraptor.yaml"
    benchmark:
        "benchmarks/varlociraptor/call/{group}.{caller}.{scatteritem}.tsv"
    shell:
        "(varlociraptor call variants {params.extra} generic --obs {params.obs}"
        " --scenario {input.scenario} {params.postprocess} {output}) 2> {log}"


rule sort_calls:
    input:
        "results/calls/{group}.{caller}.{scatteritem}.bcf",
    output:
        temp("results/calls/{group}.{caller}.{scatteritem}.sorted.bcf"),
    log:
        "logs/bcf-sort/{group}.{caller}.{scatteritem}.log",
    conda:
        "../envs/bcftools.yaml"
    resources:
        mem_mb=8000,
    shell:
        "bcftools sort --max-mem {resources.mem_mb}M --temp-dir `mktemp -d` "
        "-Ob {input} > {output} 2> {log}"


rule bcftools_concat:
    input:
        calls=get_scattered_calls(),
        indexes=get_scattered_calls(ext="bcf.csi"),
    output:
        "results/calls/{group}.{scatteritem}.bcf",
    log:
        "logs/concat-calls/{group}.{scatteritem}.log",
    params:
        extra="-a",  # TODO Check this
    wrapper:
        "v1.2.0/bio/bcftools/concat"
