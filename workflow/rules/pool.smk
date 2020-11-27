## person id of the uploaded hg19 candidates
person_id = os.path.basename(config["varvis_hg19_csv"]).split("_")[0]


rule varvis_hg19_csv_to_candidates_vcf:
    input:
        config["varvis_hg19_csv"]
    output:
        "results/candidates/hg19/{person_id}.hg19_candidates.vcf"
    conda:
        "../envs/plink.yaml"
    log:
        "results/log/{person_id}.hg19_candidates.vcf.log"
    shell:
        "bash workflow/scripts/csv2vcf.sh {input} {output} "
        "> {log} "


# rule annotatevarvis:
    # input:
        # "varvis.tsv"
        #   "varlociraptor.vcf"
    # output:
        # "results/parent_annotated/{candidates}.{poolfather}.{poolmothers}.csv"
    # shell:
        # ""
