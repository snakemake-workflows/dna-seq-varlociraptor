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
        sample=lambda wc: f'results/observations/{{group}}/{samples.loc[(samples["group"] == wc.group) & (samples["alias"] == "tumor"), "sample_name"].squeeze()}.bcf',
        contaminant=lambda wc: f'results/observations/{{group}}/{samples.loc[(samples["group"] == wc.group) & (samples["alias"] == "normal"), "sample_name"].squeeze()}.bcf',
    output:
        tsv="results/contamination/{group}.contamination_estimate.tsv",
        plot="results/contamination/{group}.contamination_estimate.json",
    log:
        "logs/varlociraptor/contamination/{group}.contamination_estimate.log",
    conda:
        "../envs/varlociraptor.yaml"
    shell:
        "varlociraptor estimate contamination "
        " --sample {input.sample} "
        " --contaminant {input.contaminant} "
        " --output-plot {output.plot} "
        " --output {output.tsv} "
        "2> {log}"
