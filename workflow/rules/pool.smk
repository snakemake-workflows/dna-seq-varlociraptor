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


## rule hg19_candidates_vcf_to_bed:
##     input:
##         "results/candidates/hg19/{person_id}.hg19_candidates.vcf"
##     output:
##         "results/candidates/hg19/{person_id}.hg19_candidates.bed"
##     shell:
##         "cat {input} "
##         "| grep -v '^#' "
##         """| awk 'BEGIN{{OFS="\t"}}{{l=length($5);print $1, $2-1, $2+l-1, $3}}' """
##         "> {output} "


rule download_hg19ToHg38_over_chain:
    output:
        "resources/hg19ToHg38.over.chain.gz" # UCSC style chromosome names
    cache:
        True
    shell:
        "curl -L http://hgdownload.soe.ucsc.edu/goldenPath/hg19/liftOver/hg19ToHg38.over.chain.gz "
        "| pigz -dc "
        "| sed 's/chr//g' " # convert chromosome names to Ensembl style
        "| pigz -9 "
        "> {output} "


rule picard_liftovervcf_hg19ToHg38:
    input:
        vcf = "results/candidates/hg19/{person_id}.hg19_candidates.vcf",
        chain = "resources/hg19ToHg38.over.chain.gz",
        fasta = "resources/genome.fasta"
    output:
        vcf="results/candidates/hg38/{person_id}.hg38_candidates.vcf",
        reject="results/candidates/hg38/{person_id}.hg38_candidates.rejected.vcf"
    log:
        "results/log/{person_id}.picard_liftovervcf.log"
    conda:
        "../envs/picard.yaml"
    shell:
        "picard -Xmx6g LiftoverVcf "
        "I={input.vcf} "
        "CHAIN={input.chain} "
        "R={input.fasta} "
        "O={output.vcf} "
        "REJECT={output.reject} "
        ">{log} 2>&1 "


rule prepocess_father_pool:
    input:
        candidates_vcf = "results/candidates/hg38/{person_id}.hg38_candidates.vcf",
        pool_bam = "results/recal/{father_pool_id}.sorted.bam",
        fasta = "resources/genome.fasta"
    output:
        "results/observations/{person_id}.{father_pool_id}.father_pool_observations.bcf"
    log:
        "results/log/{person_id}.{father_pool_id}.prepocess_father_pool.log"
    conda:
        "../envs/varlociraptor.yaml"
    shell:
        "varlociraptor preprocess variants "
        " {input.fasta} "
        "--bam {input.pool_bam} "
        "--candidates {input.candidates_vcf} "
        ">{output} "
        "2>{log} "


rule prepocess_mother_pool:
    input:
        candidates_vcf = "results/candidates/hg38/{person_id}.hg38_candidates.vcf",
        pool_bam = "results/recal/{mother_pool_id}.sorted.bam",
        fasta = "resources/genome.fasta"
    output:
        "results/observations/{person_id}.{mother_pool_id}.mother_pool_observations.bcf"
    log:
        "results/log/{person_id}.{mother_pool_id}.prepocess_mother_pool.log"
    conda:
        "../envs/varlociraptor.yaml"
    shell:
        "varlociraptor preprocess variants "
        " {input.fasta} "
        "--bam {input.pool_bam} "
        "--candidates {input.candidates_vcf} "
        ">{output} "
        "2>{log} "


rule calls:
    input:
        scenario = "config/parent_pools.scenario.yaml",
        father_pool = "results/observations/{person_id}.{father_pool_id}.father_pool_observations.bcf",
        mother_pool = "results/observations/{person_id}.{mother_pool_id}.mother_pool_observations.bcf"
    output:
        "results/calls/hg38/{person_id}.{father_pool_id}.{mother_pool_id}.calls.bcf"
    log:
        "results/log/{person_id}.{father_pool_id}.{mother_pool_id}.hg38.calls.log"
    conda:
        "../envs/varlociraptor.yaml"
    shell:
        "varlociraptor call variants "
        "generic "
        "--scenario {input.scenario} "
        "--obs "
        "father_pool={input.father_pool} "
        "mother_pool={input.mother_pool} "
        ">{output} "
        "2>{log} "


rule calls_bcf_to_vcf:
    input:
        "results/calls/hg38/{person_id}.{father_pool_id}.{mother_pool_id}.calls.bcf"
    output:
        "results/calls/hg38/{person_id}.{father_pool_id}.{mother_pool_id}.calls.vcf"
    conda:
        "../envs/bcftools.yaml"
    shell:
        "bcftools view {input} "
        ">{output} "


rule download_hg38ToHg19_over_chain:
    output:
        "resources/hg38ToHg19.over.chain.gz" # UCSC style chromosome names
    cache:
        True
    shell:
        "curl -L http://hgdownload.soe.ucsc.edu/goldenPath/hg38/liftOver/hg38ToHg19.over.chain.gz "
        "| pigz -dc "
        "| sed 's/chr//g' " # convert chromosome names to Ensembl style
        "| pigz -9 "
        "> {output} "


rule download_hg19_fasta:
    output:
        fasta = "resources/genome_hg19.fasta", # UCSC style chromosome names
        dict = "resources/genome_hg19.fasta.dict"
    cache:
        True
    conda:
        "../envs/picard.yaml"
    shell:
        "curl -L ftp://ftp.ensembl.org/pub/grch37/release-100/fasta/homo_sapiens/dna/Homo_sapiens.GRCh37.dna.primary_assembly.fa.gz "
        "| pigz -dc "
        "> {output.fasta} "
        "&& picard CreateSequenceDictionary "
        "REFERENCE={output.fasta} "
        "OUTPUT={output.dict} "


rule picard_liftovervcf_hg38ToHg19:
    input:
        vcf = "results/calls/hg38/{person_id}.{father_pool_id}.{mother_pool_id}.calls.vcf",
        chain = "resources/hg38ToHg19.over.chain.gz",
        fasta = "resources/genome_hg19.fasta"
    output:
        vcf="results/calls/hg19/{person_id}.{father_pool_id}.{mother_pool_id}.calls.vcf",
        reject="results/calls/hg19/{person_id}.{father_pool_id}.{mother_pool_id}.calls.rejected.vcf"
    log:
        "results/log/{person_id}.{father_pool_id}.{mother_pool_id}.picard_liftovervcf_hg38ToHg19.log"
    conda:
        "../envs/picard.yaml"
    shell:
        "picard -Xmx6g LiftoverVcf "
        "I={input.vcf} "
        "CHAIN={input.chain} "
        "R={input.fasta} "
        "O={output.vcf} "
        "REJECT={output.reject} "
        ">{log} 2>&1 "



# rule annotate_varvis:
#     input:
#         "varvis.tsv"
#           "varlociraptor.vcf"
#     output:
#         "results/parent_annotated/{candidates}.{father_pool_id}.{mother_pool_id}.csv"
#     shell:
#         ""

