rule render_scenario:
    input:
        local(config["calling"]["scenario"])
    output:
        report("results/scenarios/{group}.yaml", caption="../report/scenario.rst", category="Variant calling scenarios")
    log:
        "logs/render-scenario/{group}.log"
    params:
        samples = samples
    conda:
        "../envs/render_scenario.yaml"
    script:
        "../scripts/render-scenario.py"


rule varlociraptor_preprocess:
    input:
        ref="resources/genome.fasta",
        ref_idx="resources/genome.fasta.fai",
        candidates=get_candidate_calls(),
        bam="results/recal/{sample}.sorted.bam",
        bai="results/recal/{sample}.sorted.bai"
    output:
        temp("results/observations/{group}/{sample}.{caller}.{scatteritem}.bcf")
    params:
        omit_isize = "--omit-insert-size" if is_activated("primers/trimming") else ""
    log:
        "logs/varlociraptor/preprocess/{group}/{sample}.{caller}.{scatteritem}.log"
    benchmark:
        "benchmarks/varlociraptor/preprocess/{group}/{sample}.{caller}.{scatteritem}.tsv"
    conda:
        "../envs/varlociraptor.yaml"
    shell:
        "varlociraptor preprocess variants {params.omit_isize} --candidates {input.candidates} "
        "{input.ref} --bam {input.bam} --output {output} 2> {log}"


rule varlociraptor_call:
    input:
        obs=get_group_observations,
        scenario="results/scenarios/{group}.yaml"
    output:
        temp("results/calls/{group}.{caller}.{scatteritem}.bcf")
    log:
        "logs/varlociraptor/call/{group}.{caller}.{scatteritem}.log"
    params:
        obs=lambda w, input: ["{}={}".format(s, f) for s, f in zip(get_group_aliases(w), input.obs)]
    conda:
        "../envs/varlociraptor.yaml"
    benchmark:
         "benchmarks/varlociraptor/call/{group}.{caller}.{scatteritem}.tsv"
    shell:
        "varlociraptor "
        "call variants generic --obs {params.obs} "
        "--scenario {input.scenario} > {output} 2> {log}"


rule sort_calls:
    input:
       "results/calls/{group}.{caller}.{scatteritem}.bcf",
    output:
        temp("results/calls/{group}.{caller}.{scatteritem}.sorted.bcf")
    log:
        "logs/bcf-sort/{group}.{caller}.{scatteritem}.log"
    conda:
        "../envs/bcftools.yaml"
    resources:
        mem_mb=8000
    shell:
        "bcftools sort --max-mem {resources.mem_mb}M --temp-dir `mktemp -d` "
        "-Ob {input} > {output} 2> {log}"


rule bcftools_concat:
    input:
        calls = get_scattered_calls(),
        indexes = get_scattered_calls(ext=".bcf.csi"),
    output:
        "results/calls/{group}.{scatteritem}.bcf"
    log:
        "logs/condat-calls/{group}.{scatteritem}.log"
    params:
        "-a -Ob" # TODO Check this
    wrapper:
        "0.59.2/bio/bcftools/concat"
