if config["mutational_burden"]["activate"]:

    rule plot_mutational_burden:
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


    rule plot_group_vafs:
        input:
            "results/final-calls/{group}.annotated.bcf",
        output:
            report(
                "results/plots/vafs/{group}.vafs.svg",
                caption="../report/vafs.rst",
                category="Allele frequency distribution",
                subcategory="{group}",
            ),
        log:
            "logs/vafs/{group}.log",
        params:
            events=get_mutational_burden_events,
            sample=get_sample_alias,
        conda:
            "../envs/varlociraptor.yaml"
        shell:
            "(varlociraptor plot scatter "
            "--somatic-tumor-events {params.events} "
            "--sample-x {params.sample} "
            "< {input} | vl2svg > {output}) 2> {log}"