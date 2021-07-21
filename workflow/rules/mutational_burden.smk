if config["mutational_burden"]["activate"]:

    rule estimate_mutational_burden:
        input:
            "results/final-calls/{group}.annotated.bcf",
        output:
            report(
                "results/plots/mutational-burden/{group}.{sample}.{mode}.mutational-burden.svg",
                caption="../report/mutational_burden.rst",
                category="Mutational Burden",
                subcategory="{group}",
            ),
        log:
            "logs/estimate-mutational-burden/{group}.{sample}.{mode}.log",
        params:
            coding_genome_size=config["mutational_burden"]["coding_genome_size"],
            events=get_mutational_burden_events,
            sample=get_sample_alias,
        conda:
            "../envs/varlociraptor.yaml"
        shell:
            "(varlociraptor estimate mutational-burden "
            "--plot-mode {wildcards.mode} "
            "--coding-genome-size {params.coding_genome_size} "
            "--events {params.events} "
            "--sample {params.sample} "
            "< {input} | vl2svg > {output}) 2> {log}"
