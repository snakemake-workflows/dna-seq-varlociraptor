rule vembrane_table:
    input:
        bcf="results/merged-calls/{group}.{event}.fdr-controlled.bcf",
    output:
        bcf="results/tables/{group}.{event}.fdr-controlled.tsv"
    conda:
        "../envs/vembrane.yaml"
    params:
        expression="INDEX, CHROM, POS, REF, ALT[0], ANN['Consequence'], ANN['IMPACT'], ANN['SYMBOL'], ANN['Feature'], INFO['gnomad_AF']",
        prob=lambda wc: ", ".join(f"1-10**(-INFO['PROB_{x.upper()}']/10)" for x in config["calling"]["fdr-control"]["events"][wc.event]["varlociraptor"]),
        geno=lambda wc: ", ".join(f"FORMAT['AF']['{sample}']" for sample in get_group_sample_aliases(wc)),
        depth=lambda wc: ", ".join(f"FORMAT['DP']['{sample}']" for sample in get_group_sample_aliases(wc)),
        #geno="FORMAT['DP']['index'], FORMAT['AF']['index'], FORMAT['DP']['father'], FORMAT['AF']['father'], FORMAT['DP']['mother'], FORMAT['AF']['mother']",
    log:
        "logs/vembrane-table/{group}.{event}.log"
    shell:
        "vembrane table \"{params.expression}, {params.prob}, {params.geno}, {params.depth}\" {input.bcf} > {output.bcf} 2> {log}" #"1-10**(-INFO['PROB_DENOVO']/10)"


rule tsv_to_excel:
    input:
        tsv="{x}.tsv"
    output:
        xlsx="{x}.xlsx"
    conda:
        "../envs/pandas.yaml"
    script:
        "../scripts/tsv_to_xlsx.py"

rule igv_report:
    input:
        bcf="results/merged-calls/{group}.{event}.fdr-controlled.bcf",
        ref="resources/genome.fasta",
        ref_idx="resources/genome.fasta.fai",
        bams=lambda w: get_group_bams(w),
        bais=lambda w: get_group_bams(w, bai=True)
    output:
        report("results/igv-report/{group}.{event}.html", caption="../report/calls.rst", category="Variant calls")
    log:
        "logs/igv-report/{group}.{event}.log"
    params:
        os.path.join(os.path.dirname(workflow.snakefile), "resources/igv-report-template.html")
    conda:
        "../envs/igv-reports.yaml"
    shell:
        "bcftools view {input.bcf} > {output}.vcf;"
        "create_report {output}.vcf {input.ref} --flanking 100 "
        "--info-columns ANN dgiDB_drugs cosmic_LEGACY_ID --info-columns-prefixes PROB_ --sample-columns DP AF OBS"
        " --template {params} --tracks {input.bams} --output {output} --standalone 2>&1 > {log}; "
        "rm {output}.vcf"
