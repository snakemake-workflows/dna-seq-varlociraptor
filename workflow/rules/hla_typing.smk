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
        genome=genome,
        xml="results/preparation/hla.xml",
        alleles="results/preparation/hla_gen.fasta",
    output:
        temp(directory("results/orthanq-candidate-calls-temp/"))
    log:
        "logs/orthanq_candidates/candidates.log",
    threads: 8
    params:
        command="candidates",
        subcommand="hla",
        output_bcf=True 
    wrapper:
        "v7.7.0/bio/orthanq"


rule rename_candidate_bcf:
    input:
        "results/orthanq-candidate-calls-temp/"
    output:
        "results/orthanq-candidate-calls/orthanq-{locus}.hla-variants.bcf"
    log:
        "logs/scatter-candidates/rename_candidates_{locus}.log",
    shell:
        "cp {input}/{wildcards.locus}.bcf {output}"


# The scatter_candidates rule does not work without group wildcards.
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


rule orthanq_call_hla:
    input:
        haplotype_calls="results/calls/{group}.{caller}.bcf",
        haplotype_variants="results/orthanq-candidate-calls/{caller}.hla-variants.bcf",
        xml="results/preparation/hla.xml",
    output:
        directory("results/hla-typing/{group}-{caller}/{alias}")
    log:
        "logs/orthanq/{group}-{alias}-{caller}.log",
    params:
        command="call",
        subcommand="hla",
        prior="diploid",
        extra=""
    wrapper:
        "v7.7.0/bio/orthanq"
        
# TODO add other outputs (plots), fill missing inputs and commands
# TODO decide how to handle deactivation of biases in varlociraptor:
# my current favorite is a special INFO tag in the orthanq candidates.
