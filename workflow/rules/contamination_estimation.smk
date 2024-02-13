use rule bcftools_concat as bcftools_concat_all_obs_per_sample with:
    input:
        gather.calling(
            "results/observations/{{{{group}}}}/{{{{sample}}}}.{caller}.{{scatteritem}}.bcf".format(
                caller=caller
            )
        ),
    output:
        "results/observations/{group}/{sample}.bcf",
    log:
        "logs/observations/{group}/{sample}.bcftools_concat_all_obs_per_sample.log",


rule varlociraptor_estimate_contamination:
    input:
        sample="results/observations/{group}/{sample}.bcf",
        contaminant="results/observations/{group}/{sample}.bcf",
    output:
        tsv="results/contamination/{group}.{tumor_sample}.{normal_sample}_contamination_estimate.tsv",
        plot="results/contamination/{group}.{tumor_sample}.{normal_sample}_contamination_estimate.json",
    log:
        "logs/varlociraptor/contamination/{group}.{tumor_sample}.{normal_sample}_contamination_estimate.tsv",
    conda:
        "../envs/varlociraptor.yaml"
    shell:
        "varlociraptor estimate contamination "
        " --sample {input.sample} "
        " --contaminant {input.contaminant} "
        " --output-plot {output.plot} "
        " --output {output.tsv} "
        "2> {log}"
