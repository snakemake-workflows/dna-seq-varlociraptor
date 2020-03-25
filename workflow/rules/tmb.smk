if config["tmb"]["activate"]:
    rule estimate_tmb:
        input:
            "results/calls/{group}.annotated.bcf"
        output:
            "results/plots/tmb/{group}.tmb.vl.json"
        log:
            "logs/estimate-tmb/{group}.log"
        params:
            **config["tmb"]
        conda:
            "../envs/varlociraptor.yaml"
        shell:
            "varlociraptor estimate tmb "
            "--coding-genome-size {params.coding_genome_size} "
            "--somatic-tumor-events {params.somatic_events} "
            "--tumor-sample {params.tumor_sample} "
            "< {input} > {output} 2> {log}"
