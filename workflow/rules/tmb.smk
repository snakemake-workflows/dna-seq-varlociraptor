if config["tmb"]["activate"]:
    rule estimate_tmb:
        input:
            "calls/{group}.annotated.bcf"
        output:
            "plots/tmb/{group}.tmb.vl.json"
        conda:
            "../envs/varlociraptor.yaml"
        params:
            **config["tmb"]
        shell:
            "varlociraptor estimate tmb "
            "--coding-genome-size {params.coding_genome_size} "
            "--somatic-tumor-events {params.somatic_events} "
            "--tumor-sample {params.tumor_sample} "
            "< {input} > {output}"
