if config["mutational_burden"]["activate"]:

    rule calculate_covered_coding_sites:
        input:
            covered="results/regions/{group}.covered_regions.filtered.bed",
            coding="resources/coding_regions.bed.gz",
        output:
            "results/regions/{group}.covered_regions.filtered.coding.coverage_breadth.txt",
        conda:
            "../envs/awk_bedtools.yaml"
        log:
            "logs/regions/{group}.covered_regions.filtered.coding.coverage_breadth.log",
        shell:
            "( bedtools intersect -a {input.covered} -b {input.coding} | "
            "  awk '{{ covered += $3 - $2 }} END {{ print covered }}' >{output} "
            ") 2> {log} "

    rule estimate_mutational_burden:
        input:
            calls="results/final-calls/{group}.annotated.bcf",
            coverage_breadth="results/regions/{group}.covered_regions.filtered.coding.coverage_breadth.txt",
        output:
            report(
                "results/plots/mutational-burden/{group}.{sample}.{mode}.mutational-burden.svg",
                caption="../report/mutational_burden.rst",
                category="Mutational Burden",
                subcategory="{group}",
                labels={"sample": "{sample}", "mode": "{mode}"},
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
