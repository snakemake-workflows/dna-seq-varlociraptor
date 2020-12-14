rule map_reads:
    input:
        reads=get_map_reads_input,
        idx=rules.bwa_index.output
    output:
        temp("results/mapped/{sample}.sorted.bam")
    log:
        "logs/bwa_mem/{sample}.log"
    params:
        index=lambda w, input: os.path.splitext(input.idx[0])[0],
        extra=get_read_group,
        sort="samtools",
        sort_order="coordinate"
    threads: 8
    wrapper:
        "0.56.0/bio/bwa/mem"


rule mark_duplicates:
    input:
        "results/mapped/{sample}.sorted.bam"
    output:
        bam=temp("results/dedup/{sample}.sorted.bam"),
        metrics="results/qc/dedup/{sample}.metrics.txt"
    log:
        "logs/picard/dedup/{sample}.log"
    params:
        extra = config["mark_duplicates"].get("params", ""),
    wrapper:
        "master/bio/picard/markduplicateswithmatecigar"


############
### BQSR ###
############

rule indel_realign_create:
    input:
        unpack(lambda wildcards: get_bqsr_input(wildcards, "indel_realign")),
        ref = expand(rules.get_genome.output, **config['ref']),
        known = expand(rules.remove_iupac_codes.output, **config['ref']),
        tbi = expand(str(rules.remove_iupac_codes.output) + ".tbi", **config['ref']),
    output:
        intervals = temp("temp/indel_realign/create/{sample}.intervals")
    log:
        "logs/indel_realign/create/{sample}.log"
    params:
        extra = config["indel_realignment"].get("params", {}).get("target_creator",""),
    threads: 16
    resources:
        mem_mb = 10 * 1024,
        time_min = 10 * 60,
    wrapper:
        "master/bio/gatk3/realignertargetcreator"

ruleorder: indel_realign > bam_index
rule indel_realign:
    input:
        unpack(lambda wildcards: get_bqsr_input(wildcards, "indel_realign")),
        ref = expand(rules.get_genome.output, **config['ref']),
        known = expand(rules.remove_iupac_codes.output, **config['ref']),
        tbi = expand(str(rules.remove_iupac_codes.output) + ".tbi", **config['ref']),
        target_intervals = rules.indel_realign_create.output,
    output:
        bam = temp("temp/indel_realign/{sample}.bam"),
        bai = temp("temp/indel_realign/{sample}.bai"),
    log:
        "logs/indel_realign/{sample}.log"
    params:
        extra = config["indel_realignment"].get("params", {}).get("indel_realigner",""),
    threads: 16
    resources:
        mem_mb = 10 * 1024,
        time_min = 10 * 60,
    wrapper:
        "master/bio/gatk3/indelrealigner"


rule calmd:
    input:
        aln = lambda wildcards: get_bqsr_input(wildcards, "calmd")["bam"],
        ref = expand(rules.get_genome.output, **config['ref']),
    output:
        bam = temp("temp/bqsr/calmd/{sample}.bam"),
    log:
        "logs/bqsr/calmd/{sample}.log"
    params:
        extra = "",
    wrapper:
        "master/bio/samtools/calmd"


rule create_gatk_bqsr:
    input:
        unpack(lambda wildcards: get_bqsr_input(wildcards, "gatk_bqsr")),
        ref=expand(rules.get_genome.output, **config['ref']),
        ref_dict=expand(rules.genome_dict.output, **config['ref']),
        ref_fai=expand(rules.genome_faidx.output, **config['ref']),
        known=expand(rules.remove_iupac_codes.output, **config['ref']),
        tbi=expand(str(rules.remove_iupac_codes.output) + ".tbi", **config['ref']),
    output:
        recal_table=temp("temp/bqsr/gatk_bqsr/{sample}.grp"),
    params:
        extra=config["bqsr"]["gatk_bqsr"].get("params", {}).get("create_bqsr", ""),
    log:
        "logs/bqsr/gatk_bqsr/create/{sample}.log"
    threads: 8
    resources:
        mem_mb = 10 * 1024,
        time_min = 10 * 60,
    wrapper:
        "0.67.0/bio/gatk/baserecalibratorspark"


ruleorder: apply_gatk_bqsr > bam_index
rule apply_gatk_bqsr:
    input:
        unpack(lambda wildcards: get_bqsr_input(wildcards, "gatk_bqsr")),
        ref=expand(rules.get_genome.output, **config['ref']),
        ref_dict=expand(rules.genome_dict.output, **config['ref']),
        ref_fai=expand(rules.genome_faidx.output, **config['ref']),
        recal_table=rules.create_gatk_bqsr.output.recal_table,
    output:
        bam = temp("temp/bqsr/gatk_bqsr/{sample}.bam"),
        bai = temp("temp/bqsr/gatk_bqsr/{sample}.bai"),
    log:
        "logs/bqsr/gatk/gatk_bqsr/apply/{sample}.log"
    params:
        extra = config["bqsr"]["gatk_bqsr"].get("params", {}).get("apply_bqsr", ""),
    threads: 1
    resources:
        mem_mb = 10 * 1024,
        time_min = 12 * 60,
    wrapper:
        "0.67.0/bio/gatk/applybqsr"


rule cram:
    input:
        unpack(lambda wildcards: get_bqsr_input(wildcards, "cram")),
        ref = expand(rules.get_genome.output, **config['ref']),
    output:
        cram = protected("results/cram/{sample}.cram"),
    log:
        "logs/cram/{sample}.log"
    params:
        lambda wildcards, input: f"--reference {input.ref} --output-fmt CRAM"
    wrapper:
        "0.67.0/bio/samtools/view"
