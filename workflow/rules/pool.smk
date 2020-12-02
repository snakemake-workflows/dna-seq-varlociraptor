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


## ##FORMAT=<ID=OBS,Number=A,Type=String,Description="Summary of observations. Each entry is encoded as CBTSO,
## with C being a count,
## B being the posterior odds for the alt allele (see below),
## T being the type of alignment, encoded as s=single end and p=paired end,
## S being the strand that supports the observation (+, -, or * for both), and
## O being the read orientation (> = F1R2, < = F2R1, * = unknown, ! = non standard, e.g. R1F2).
## Posterior odds (B) for alt allele of each fragment are given as extended Kass Raftery scores:
##      N=none, E=equal, B=barely, P=positive, S=strong, V=very strong (lower case if probability for correct mapping
##      of fragment is <95%). Thereby we extend Kass Raftery scores
##      with a term for equality between the evidence of the two alleles (E=equal).">
##
## example:
## 1Np+<28Np-<22Np+>15Np->4Ns+*4Ns-*
## 1Np+< 28Np-< 22Np+> 15Np-> 4Ns+* 4Ns-*
##
## 1 mal Referenzbase mit forward readpair
## 28 mal Referenzbase mit reverse readpair
## 22 mal Referenzbase mit forward readpair
## 15 mal Referenzbase mit reverse readpair
## 4 forward Reads haben hier einen split
## 4 reverse Reads haben hier einen split


rule vembrane:
    input:
        "results/calls/hg38/{person_id}.{father_pool_id}.{mother_pool_id}.calls.vcf"
    output:
        "results/calls/hg38/AF/{person_id}.{father_pool_id}.{mother_pool_id}.AF.tsv"
    conda:
        "../envs/vembrane.yaml"
    shell:
        "vembrane table "
        "--header 'ID, FATHER_N_REF, FATHER_N_ALT, FATHER_PROB, MOTHER_N_REF, MOTHER_N_ALT, MOTHER_PROB' "
        """"ID, """
        """sum(map(int, re.findall('(\d+)[N]', FORMAT['OBS']['father_pool']))), """ 
        """sum(map(int, re.findall('(\d+)[VS]', FORMAT['OBS']['father_pool']))), """
        ## """1.0 * sum(map(int, re.findall('(\d+)[VS]', FORMAT['OBS']['father_pool']))) / sum(map(int, re.findall('(\d+)[NVS]', FORMAT['OBS']['father_pool']))), """
        """10 ** (-INFO['PROB_FATHER_ONLY']/10), """
        """sum(map(int, re.findall('(\d+)[N]', FORMAT['OBS']['mother_pool']))), """
        """sum(map(int, re.findall('(\d+)[VS]', FORMAT['OBS']['mother_pool']))), """
        ## """1.0 * sum(map(int, re.findall('(\d+)[VS]', FORMAT['OBS']['mother_pool']))) / sum(map(int, re.findall('(\d+)[NVS]', FORMAT['OBS']['mother_pool']))), """
        """10 ** (-INFO['PROB_MOTHER_ONLY']/10)" """
        "{input} "
        "| dos2unix "
        ">{output}.temp "
        "&& ("
        """echo -e "$(head -1 {output}.temp)\tFATHER_AF\tMOTHER_AF" """
        "&& tail -n+2 {output}.temp "
        """| awk 'BEGIN{{OFS="\t"}}{{if ($2==0) {{FATHER_AF="NA"}} else {{FATHER_AF=$3/($2+$3)}}; if ($5==0) {{MOTHER_AF="NA"}} else {{MOTHER_AF=$6/($5+$6)}}; print $0, FATHER_AF, MOTHER_AF}}'"""
        ") >{output} "
        "&& rm {output}.temp "


# rule annotate_varvis:
#     input:
#         "varvis.tsv"
#           "varlociraptor.vcf"
#     output:
#         "results/parent_annotated/{candidates}.{father_pool_id}.{mother_pool_id}.csv"
#     shell:
#         ""

