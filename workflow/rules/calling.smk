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
        candidates="results/candidate-calls/{group}.{caller}.bcf",
        bam="results/recal/{sample}.sorted.bam",
        bai="results/recal/{sample}.sorted.bam.bai"
    output:
        temp("results/observations/{group}/{sample}.{caller}.bcf")
    params:
        omit_isize = "--omit-insert-size" if is_activated("primers/trimming") else ""
    log:
        "logs/varlociraptor/preprocess/{group}/{sample}.{caller}.log"
    threads: workflow.cores
    benchmark:
        "benchmarks/varlociraptor/preprocess/{group}/{sample}.{caller}.tsv"
    conda:
        "../envs/varlociraptor.yaml"
    shell:
        "varlociraptor preprocess variants {params.omit_isize} --candidates {input.candidates} "
        "{input.ref} --bam {input.bam} --output {output} --threads {threads} 2> {log}"


rule sort_observations:
    input:
       "results/observations/{group}/{sample}.{caller}.bcf"
    output:
        "results/observations/{group}/{sample}.{caller}.sorted.bcf"
    log:
        "logs/sort-observations/{group}.{sample}.{caller}.log"
    conda:
        "../envs/bcftools.yaml"
    shell:
        "bcftools sort --temp-dir $TMPDIR "
        "-Ob {input} > {output} 2> {log}"


rule varlociraptor_call:
    input:
        obs=get_group_observations,
        scenario="results/scenarios/{group}.yaml"
    output:
        temp("results/calls/{group}.{caller}.bcf")
    log:
        "logs/varlociraptor/call/{group}.{caller}.log"
    params:
        obs=lambda w, input: ["{}={}".format(s, f) for s, f in zip(get_group_aliases(w), input.obs)]
    threads: workflow.cores
    conda:
        "../envs/varlociraptor.yaml"
    shell:
        "varlociraptor "
        "call variants --threads {threads} generic --obs {params.obs} "
        "--scenario {input.scenario} > {output} 2> {log}"


rule sort_calls:
    input:
       "results/calls/{group}.{caller}.bcf",
    output:
        temp("results/calls/{group}.{caller}.sorted.bcf")
    log:
        "logs/bcf-sort/{group}.{caller}.log"
    conda:
        "../envs/bcftools.yaml"
    shell:
        "bcftools sort --temp-dir $TMPDIR "
        "-Ob {input} > {output} 2> {log}"


rule bcftools_concat:
    input:
        calls = expand(
            "results/calls/{{group}}.{caller}.sorted.bcf",
            caller=caller
        ),
        indexes = expand(
            "results/calls/{{group}}.{caller}.sorted.bcf.csi",
            caller=caller
        )
    output:
        "results/calls/{group}.bcf"
    log:
        "logs/condat-calls/{group}.log"
    params:
        "-a -Ob" # TODO Check this
    wrapper:
        "0.59.2/bio/bcftools/concat"
