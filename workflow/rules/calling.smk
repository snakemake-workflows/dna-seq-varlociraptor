rule render_scenario:
    input:
        config["calling"]["scenario"]
    output:
        report("scenarios/{group}.yaml", caption="../report/scenario.rst", category="Variant calling scenarios")
    params:
        samples = samples
    conda:
        "../envs/render_scenario.yaml"
    script:
        "../scripts/render-scenario.py"

rule varlociraptor_preprocess:
    input:
        ref="refs/genome.fasta",
        ref_idx="refs/genome.fasta.fai",
        candidates="candidate-calls/{group}.{caller}.bcf",
        bam="recal/{sample}.sorted.bam",
        bai="recal/{sample}.sorted.bam.bai"
    output:
        "observations/{group}/{sample}.{caller}.bcf"
    log:
        "logs/varlociraptor/preprocess/{group}/{sample}.{caller}.log"
    conda:
        "../envs/varlociraptor.yaml"
    shell:
        "varlociraptor preprocess variants --candidates {input.candidates} "
        "{input.ref} --bam {input.bam} --output {output} 2> {log}"

rule varlociraptor_call:
    input:
        obs=get_group_observations,
        scenario="scenarios/{group}.yaml"
    output:
        temp("calls/{group}.{caller}.bcf")
    log:
        "logs/varlociraptor/call/{group}.{caller}.log"
    params:
        obs=lambda w, input: ["{}={}".format(s, f) for s, f in zip(get_group_aliases(w), input.obs)]
    conda:
        "../envs/varlociraptor.yaml"
    shell:
        "varlociraptor "
        "call variants generic --obs {params.obs} "
        "--scenario {input.scenario} > {output} 2> {log}"

rule bcftools_concat:
    input:
        calls = expand(
            "calls/{{group}}.{caller}.bcf",
            caller=caller
        ),
        indexes = expand(
            "calls/{{group}}.{caller}.bcf.csi",
            caller=caller
        )
    output:
        "calls/{group}.bcf"
    params:
        "-a -Ob" # Check this
    wrapper:
        "0.36.0/bio/bcftools/concat"
