# dowloads latest IMGT/HLA database version
rule get_hla_genes_and_xml:
    output:
        genes="results/preparation/hla_gen.fasta",
        xml="results/preparation/hla.xml.zip",
    log:
        "logs/get_hla_genes_and_xml.log",
    params:
        genes_link="ftp://ftp.ebi.ac.uk/pub/databases/ipd/imgt/hla/hla_gen.fasta",
        xml_link="ftp://ftp.ebi.ac.uk/pub/databases/ipd/imgt/hla/xml/hla.xml.zip",
    shell:
        "wget -c {params.genes_link} -O {output.genes} && "
        "wget -c {params.xml_link} -O {output.xml} 2> {log}"


rule unzip_xml:
    input:
        "results/preparation/hla.xml.zip",
    output:
        xml="results/preparation/hla.xml",
    log:
        "logs/unzip_xml.log",
    params:
        path_to_unzip=lambda wc, output: os.path.dirname(output.xml),
    shell:
        "unzip -o {input} -d {params.path_to_unzip}"


rule orthanq_candidate_variants:
    input:
        hla_genes="results/preparation/hla_gen.fasta",
        xml="results/preparation/hla.xml",
        genome=genome,
    output:
        bcfs=expand(
            "results/orthanq-candidate-calls/{caller}.hla-variants.bcf",
            caller=[f"orthanq-{locus}" for locus in config["hla_typing"].get("loci")],
        ),
    log:
        "logs/orthanq-candidates.log",
    conda:
        "../envs/orthanq.yaml"
    params:
        output_folder=lambda wc, output: os.path.dirname(output.bcfs[0]),
    threads: 8
    shell:
        "orthanq candidates hla --alleles {input.hla_genes} --genome {input.genome} --xml {input.xml} "
        "--threads {threads} --output-bcf --output {params.output_folder} 2> {log}"


# since the scatter_candidates rule does not work without group wildcards.
rule scatter_candidates_orthanq:
    input:
        "results/orthanq-candidate-calls/{orthanq_locus}.hla-variants.bcf",
    output:
        scatter.calling(
            "results/orthanq-candidate-calls/{{orthanq_locus}}.hla-variants.{scatteritem}.bcf"
        ),
    log:
        "logs/scatter-candidates/{orthanq_locus}.log",
    conda:
        "../envs/rbt.yaml"
    shell:
        "rbt vcf-split {input} {output}"


rule gather_annotated_calls_orthanq:
    input:
        calls=gather.calling("results/calls/{{group}}.{{caller}}.{scatteritem}.bcf"),
        idx=gather.calling("results/calls/{{group}}.{{caller}}.{scatteritem}.bcf.csi"),
    output:
        "results/calls/{group}.{caller}.bcf",
    log:
        "logs/gather-hla-variants/{group}.{caller}.log",
    params:
        extra="-a",
    group:
        "annotation"
    wrapper:
        "v2.3.2/bio/bcftools/concat"


rule orthanq_call:
    input:
        haplotype_variants="results/orthanq-candidate-calls/{caller}.hla-variants.bcf",
        haplotype_calls="results/calls/{group}.{caller}.bcf",
        xml="results/preparation/hla.xml",
    output:  #orthanq uses underscore for a separator of sample/group name and locus name.
        table="results/hla-typing/{group}-{caller}/{alias}/{alias}-{caller}.csv",
        three_field_solutions="results/hla-typing/{group}-{caller}/{alias}/3_field_solutions.json",
        two_field_solutions="results/hla-typing/{group}-{caller}/{alias}/2_field_solutions.json",
        final_solution="results/hla-typing/{group}-{caller}/{alias}/final_solution.json",
        lp_solution="results/hla-typing/{group}-{caller}/{alias}/lp_solution.json",
        two_field_table="results/hla-typing/{group}-{caller}/{alias}/2-field.csv",
        g_groups="results/hla-typing/{group}-{caller}/{alias}/G_groups.csv",
    log:
        "logs/orthanq/{group}-{alias}-{caller}.log",
    conda:
        "../envs/orthanq.yaml"
    shell:
        "orthanq call hla --haplotype-variants {input.haplotype_variants} --xml {input.xml} --sample {wildcards.alias} "
        " --haplotype-calls {input.haplotype_calls} --prior diploid --output {output.table} 2> {log}"


# TODO add other outputs (plots), fill missing inputs and commands
# TODO decide how to handle deactivation of biases in varlociraptor:
# my current favorite is a special INFO tag in the orthanq candidates.
