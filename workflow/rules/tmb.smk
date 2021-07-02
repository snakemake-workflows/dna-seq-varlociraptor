if config["tmb"]["activate"]:

    rule estimate_tmb:
        input:
            "results/final-calls/{group}.annotated.bcf",
        output:
            "results/plots/tmb/{group}.{sample}.{mode}.tmb.vl.json",
        log:
            "logs/estimate-tmb/{group}.{sample}.{mode}.log",
        params:
            coding_genome_size=config["tmb"]["coding_genome_size"],
            sample=get_tmb_params,
        conda:
            "../envs/varlociraptor.yaml"
        shell:
            "varlociraptor estimate tmb "
            "--plot-mode {wildcards.mode} "
            "--coding-genome-size {params.coding_genome_size} "
            "--somatic-tumor-events {params.sample[somatic_events]} "
            "--tumor-sample {params.sample[tumor_sample]} "
            "< {input} > {output} 2> {log}"
