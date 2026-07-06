if config["mutational_burden"]["activate"]:

    rule calculate_covered_coding_sites:
        input:
            covered="results/regions/{group}.covered_regions.filtered.bed",
            coding="resources/coding_regions.bed.gz",
        output:
            "results/regions/{group}.covered_regions.filtered.coding.coverage_breadth.txt",
        log:
            "logs/regions/{group}.covered_regions.filtered.coding.coverage_breadth.log",
        conda:
            "../envs/awk_bedtools.yaml"
        shell:
            "( bedtools intersect -a {input.covered} -b {input.coding} | "
            "  awk '{{ covered += $3 - $2 }} END {{ print covered }}' >{output} "
            ") 2> {log} "

    rule estimate_mutational_burden:
        input:
            calls="results/final-calls/{group}/{group}.variants.annotated.bcf",
            coverage_breadth="results/regions/{group}.covered_regions.filtered.coding.coverage_breadth.txt",
        output:
            report(
                "results/plots/mutational-burden/{group}.{alias}.{mode}.mutational-burden.svg",
                caption="../report/mutational_burden.rst",
                category="Mutational Burden",
                subcategory="{group}",
                labels={"alias": "{alias}", "mode": "{mode}"},
            ),
        log:
            "logs/estimate-mutational-burden/{group}.{alias}.{mode}.log",
        conda:
            "../envs/varlociraptor.yaml"
        params:
            events=" ".join(config["mutational_burden"]["events"]),
        shell:
            "(varlociraptor estimate mutational-burden "
            "--mode {wildcards.mode} "
            "--coding-genome-size $( cat {input.coverage_breadth} ) "
            "--events {params.events} "
            "--sample {wildcards.alias} "
            "< {input.calls} | vl2svg > {output}) 2> {log}"
