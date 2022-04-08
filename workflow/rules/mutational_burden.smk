if config["mutational_burden"]["activate"]:

    rule calculate_covered_genomic_sites:
        input:
            "results/regions/{group}.covered_regions.filtered.bed",
        output:
            "results/regions/{group}.covered_regions.filtered.coverage_breadth.txt",
        log:
            "logs/regions/{group}.covered_regions.filtered.coverage_breadth.log",
        shell:
            "awk '{ covered += $3 - $2 } END { print covered }' {input} >{output} 2> {log} "

    rule estimate_mutational_burden:
        input:
            calls="results/final-calls/{group}.annotated.bcf",
            coverage_breadth="results/regions/{group}.covered_regions.filtered.coverage_breadth.txt",
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
            events=get_mutational_burden_events,
            sample=get_sample_alias,
        conda:
            "../envs/varlociraptor.yaml"
        shell:
            "(varlociraptor estimate mutational-burden "
            "--plot-mode {wildcards.mode} "
            "--coding-genome-size $( cat {input.coverage_breadth} ) "
            "--events {params.events} "
            "--sample {params.sample} "
            "< {input.calls} | vl2svg > {output}) 2> {log}"
