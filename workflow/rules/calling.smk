import io


rule render_scenario:
    input:
        template=local(config["calling"]["scenario"]),
    output:
        report(
            "results/scenarios/{group}.yaml",
            caption="../report/scenario.rst",
            category="Variant calling scenarios",
            labels={"sample group": "{group}"},
        ),
    log:
        "logs/render-scenario/{group}.log",
    params:
        samples=lambda wc: samples[samples["group"] == wc.group],
        annotation=lambda wc: group_annotation.loc[wc.group],
    conda:
        None
    template_engine:
        "yte"


rule varlociraptor_alignment_properties:
    input:
        ref=genome,
        ref_idx=genome_fai,
        bam=get_sample_bam,
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
        ref=genome,
        ref_idx=genome_fai,
        candidates=get_candidate_calls(),
        bam=get_sample_bam,
        bai=lambda wc: get_sample_bam(wc, bai=True),
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
        temp("results/calls/{group}.{caller}.{scatteritem}.unsorted.bcf"),
    log:
        "logs/varlociraptor/call/{group}.{caller}.{scatteritem}.log",
    params:
        obs=lambda w, input: [
            "{}={}".format(s, f) for s, f in zip(get_group_aliases(w.group), input.obs)
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
        "results/calls/{group}.{caller}.{scatteritem}.unsorted.bcf",
    output:
        temp("results/calls/{group}.{caller}.{scatteritem}.bcf"),
    params:
        # Set to True, in case you want uncompressed BCF output
        uncompressed_bcf=False,
        # Extra arguments
        extras="",
    log:
        "logs/bcf-sort/{group}.{caller}.{scatteritem}.log",
    resources:
        mem_mb=8000,
    wrapper:
        "v1.21.0/bio/bcftools/sort"


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
        "v1.21.0/bio/bcftools/concat"
