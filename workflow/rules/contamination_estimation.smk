use rule bcftools_concat as bcftools_concat_all_obs_per_sample with:
    input:
        calls=get_scattered_obs(),
        indexes=get_scattered_obs(ext="bcf.csi"),
    output:
        "results/observations/{group}/{sample}.freebayes.bcf",
    log:
        "logs/observations/{group}/{sample}.bcftools_concat_all_obs_per_sample.log",


rule varlociraptor_estimate_contamination:
    input:
        sample=lambda wc: f'results/observations/{{group}}/{samples.loc[(samples["group"] == wc.group) & (samples["alias"] == "tumor"), "sample_name"].squeeze()}.freebayes.bcf',
        contaminant=lambda wc: f'results/observations/{{group}}/{samples.loc[(samples["group"] == wc.group) & (samples["alias"] == "normal"), "sample_name"].squeeze()}.freebayes.bcf',
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
