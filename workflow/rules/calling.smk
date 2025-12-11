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
        bam="results/recal/{sample}.bam",
        bai="results/recal/{sample}.bai",
    output:
        "results/alignment-properties/{group}/{sample}.json",
    log:
        "logs/varlociraptor/estimate-alignment-properties/{group}/{sample}.log",
    conda:
        "../envs/varlociraptor.yaml"
    shell:
        "varlociraptor estimate alignment-properties {input.ref} --bams {input.bam} > {output} 2> {log}"


rule varlociraptor_preprocess:
    input:
        ref=genome,
        ref_idx=genome_fai,
        candidates=access.multi(get_candidate_calls),
        bam="results/recal/{sample}.bam",
        bai="results/recal/{sample}.bai",
        alignment_props=get_alignment_props,
    output:
        "results/observations/{caller}/{group}/{sample}.{scatteritem}.bcf",
    params:
        extra=lambda wc: get_varlociraptor_params(
            wc, config["params"]["varlociraptor"]["preprocess"]
        ),
    log:
        "logs/varlociraptor/preprocess/{group}/{sample}.{caller}.{scatteritem}.log",
    benchmark:
        "benchmarks/varlociraptor/preprocess/{group}/{sample}.{caller}.{scatteritem}.tsv"
    conda:
        "../envs/varlociraptor.yaml"
    resources:
        mem_mb=lambda wc, input, attempt: attempt * 3500,  # usually, this should stay below 3500, but we have seen RNAseq cases where this goes up to 6.5 GB
    shell:
        "varlociraptor preprocess variants --candidates {input.candidates} {params.extra} "
        "--alignment-properties {input.alignment_props} {input.ref} --bam {input.bam} --output {output} "
        "2> {log}"


rule varlociraptor_call:
    input:
        obs=get_group_observations,
        scenario="results/scenarios/{group}.yaml",
    output:
        temp("results/calls/varlociraptor/{caller}/{group}/{group}.{scatteritem}.unsorted.bcf"),
    log:
        "logs/varlociraptor/call/{caller}/{group}.{scatteritem}.log",
    params:
        obs=get_varlociraptor_obs_args,
        # Add -o as workaround to separate info field entries from subcommand
        extra=lambda wc: get_varlociraptor_params(
            wc, config["params"]["varlociraptor"]["call"]
        )
        + " -o /dev/stdout",
        postprocess=(
            ">"
            if not config["calling"].get("infer_genotypes")
            else "| varlociraptor genotype >"
        ),
    conda:
        "../envs/varlociraptor.yaml"
    benchmark:
        "benchmarks/varlociraptor/call/{group}/{group}.{caller}.{scatteritem}.tsv"
    shell:
        "(varlociraptor call variants {params.extra} generic --obs {params.obs}"
        " --scenario {input.scenario} {params.postprocess} {output}) 2> {log}"


rule sort_calls:
    input:
        "results/calls/varlociraptor/{caller}/{group}/{group}.{scatteritem}.unsorted.bcf",
    output:
        temp("results/calls/varlociraptor/{caller}/{group}/{group}.{scatteritem}.bcf"),
    params:
        # Set to True, in case you want uncompressed BCF output
        uncompressed_bcf=False,
        # Extra arguments
        extras="",
    log:
        "logs/bcf-sort/{group}/{group}.{caller}.{scatteritem}.log",
    resources:
        mem_mb=8000,
    wrapper:
        "v2.6.0/bio/bcftools/sort"


rule bcftools_concat:
    input:
        calls=get_scattered_calls(),
        indexes=get_scattered_calls(ext="bcf.csi"),
    output:
        "results/calls/varlociraptor/{group}/{group}.{calling_type}.{scatteritem}.bcf",
    log:
        "logs/concat-calls/{group}/{group}.{calling_type}.{scatteritem}.log",
    params:
        extra="-a",  # TODO Check this
    wrapper:
        "v2.3.2/bio/bcftools/concat"
