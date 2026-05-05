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
        "wget -c {params.genes_link} -O {output.genes} 2> {log} && "
        "wget -c {params.xml_link} -O {output.xml} 2>> {log}"

rule unzip_xml:
    input:
        "results/preparation/hla.xml.zip",
    output:
        xml="results/preparation/hla.xml",
    log:
        "logs/unzip_xml.log",
    conda:
        "../envs/unzip.yaml"
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
        directory("results/orthanq-candidate-calls-temp")
    log:
        "logs/orthanq_candidates/candidates.log",
    conda:
        "../envs/orthanq-deps.yaml" #only required for testing purposes
    threads: 64
    params:
        command="candidates",
        subcommand="hla",
        output_bcf=True 
    # wrapper: #commented out for testing purposes
    #     "v7.7.0/bio/orthanq"
    shell:
        "../../orthanq/target/release/orthanq candidates hla --alleles {input.alleles} --genome {input.genome} --xml {input.xml} --threads {threads} --output-bcf --output {output} 2> {log}"


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
        # output_dir=directory("results/hla-typing/{group}-{caller}/{alias}/{prior}"),
        predictions="results/hla-typing/{group}-{caller}/{alias}/{prior}/predictions.csv",
        G_groups="results/hla-typing/{group}-{caller}/{alias}/{prior}/G_groups.csv",
        three_field="results/hla-typing/{group}-{caller}/{alias}/{prior}/3_field_solutions.json",
        two_field="results/hla-typing/{group}-{caller}/{alias}/{prior}/2_field_solutions.json",
        best_solution="results/hla-typing/{group}-{caller}/{alias}/{prior}/best_solution.json",
        arrow_plot="results/hla-typing/{group}-{caller}/{alias}/{prior}/arrow_plot.json"
    log:
        "logs/orthanq/{group}-{alias}-{caller}-{prior}.log",
    # wildcard_constraints:
    #     alias="[^/]+",
    #     prior="[^/]+",
    #     caller="[^/]+",
    #     group="[^/]+"
    params:
        command="call",
        subcommand="hla",
        # prior=lambda wc: lambda wc: "diploid" if wc.alias == "normal" else "diploid-subclonal", # todo: activate when diploid-subclonal is available.
        prior= "diploid",
        extra="",
        outdir=lambda wc, output: os.path.dirname(output.three_field),
        events=lambda wc: get_events,
        allele_freqs="resources/allele_frequencies.csv"
    # wrapper:
    #     "v7.7.0/bio/orthanq"
    conda:
        "../envs/vega.yaml"
    shell:
        """
        ../../orthanq/target/release/orthanq call hla \
         --haplotype-variants {input.haplotype_variants} \
         --xml {input.xml} --haplotype-calls {input.haplotype_calls} \
         --prior {params.prior} --output {params.outdir} \
         --sample {wildcards.alias} \
         --allele-freqs {params.allele_freqs} \
         --fast \
         --events {params.events}  2> {log}
        """

rule convert_orthanq_plots:
    input:
        three_field="results/hla-typing/{group}-{caller}/{alias}/{prior}/3_field_solutions.json",
        two_field="results/hla-typing/{group}-{caller}/{alias}/{prior}/2_field_solutions.json",
        best_solution="results/hla-typing/{group}-{caller}/{alias}/{prior}/best_solution.json",
        arrow_plot="results/hla-typing/{group}-{caller}/{alias}/{prior}/arrow_plot.json"
    output:
        three_field_html="results/hla-typing/{group}-{caller}/{alias}/{prior}/3_field_solutions.html",
        three_field_svg="results/hla-typing/{group}-{caller}/{alias}/{prior}/3_field_solutions.svg",
        two_field_html="results/hla-typing/{group}-{caller}/{alias}/{prior}/2_field_solutions.html",
        two_field_svg="results/hla-typing/{group}-{caller}/{alias}/{prior}/2_field_solutions.svg",
        arrow_plot_html="results/hla-typing/{group}-{caller}/{alias}/{prior}/arrow_plot.html",
        arrow_plot_svg="results/hla-typing/{group}-{caller}/{alias}/{prior}/arrow_plot.svg",
        best_solution_html="results/hla-typing/{group}-{caller}/{alias}/{prior}/best_solution.html",
        best_solution_svg="results/hla-typing/{group}-{caller}/{alias}/{prior}/best_solution.svg",
    log:
        "logs/convert_orthanq_plots/{group}-{alias}-{caller}-{prior}.log",
    conda:
        "../envs/vega.yaml"
    shell:
        "vl-convert vl2html --input {input.three_field} --output {output.three_field_html} 2> {log} && "
        "vl2svg {input.three_field} {output.three_field_svg} 2> {log} && "
        "vl-convert vl2html --input {input.two_field} --output {output.two_field_html} 2> {log} && "
        "vl2svg {input.two_field} {output.two_field_svg} 2> {log} && "
        "vl-convert vl2html --input {input.best_solution} --output {output.best_solution_html} 2> {log} && "
        "vl2svg {input.best_solution} {output.best_solution_svg} 2> {log} && "
        "vl-convert vl2html --input {input.arrow_plot} --output {output.arrow_plot_html} 2> {log} && "
        "vl2svg {input.arrow_plot} {output.arrow_plot_svg} 2> {log}"


# TODO: decide how to handle deactivation of biases in varlociraptor, decide if --omit-mapq-adjustment is required.
# TODO: use wrappers (needs to be configured to be compliant with the latest orthanq outputs) for orthanq, make fast mode optional (remove orthanq-deps.yaml)
# TODO: enable subclonal typing with Orthanq using diploid-subclonal prior. 
