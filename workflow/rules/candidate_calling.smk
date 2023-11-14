rule freebayes:
    input:
        ref=genome,
        ref_idx=genome_fai,
        regions="results/regions/{group}.expanded_regions.filtered.bed",
        # you can have a list of samples here
        alns=lambda w: get_group_bams(w),
        idxs=lambda w: get_group_bams(w, bai=True),
    output:
        "results/candidate-calls/{group}.freebayes.bcf",
    log:
        "logs/freebayes/{group}.log",
    params:
        # genotyping is performed by varlociraptor, hence we deactivate it in freebayes by 
        # always setting --pooled-continuous
        extra="--pooled-continuous --min-alternate-count {} --min-alternate-fraction {}".format(
            1 if is_activated("calc_consensus_reads") else 2,
            config["params"]["freebayes"].get("min_alternate_fraction", "0.05"),
        ),
    threads: max(workflow.cores - 1, 1)  # use all available cores -1 (because of the pipe) for calling
    wrapper:
        "v2.7.0/bio/freebayes"


rule delly:
    input:
        ref=genome,
        ref_idx=genome_fai,
        alns=lambda w: get_group_bams(w),
        index=lambda w: get_group_bams(w, bai=True),
        exclude=get_delly_excluded_regions(),
    output:
        "results/candidate-calls/{group}.delly.bcf",
    log:
        "logs/delly/{group}.log",
    params:
        extra=config["params"].get("delly", ""),
    threads: lambda _, input: len(input.alns)  # delly parallelizes over the number of samples
    wrapper:
        "v2.3.2/bio/delly"


# Delly breakends lead to invalid BCFs after VEP annotation (invalid RLEN). Therefore we exclude them for now.
rule fix_delly_calls:
    input:
        "results/candidate-calls/{group}.delly.bcf",
    output:
        "results/candidate-calls/{group}.delly.no_bnds.bcf",
    log:
        "logs/fix_delly_calls/{group}.log",
    conda:
        "../envs/bcftools.yaml"
    shell:
        """bcftools view -e 'INFO/SVTYPE="BND"' {input} -Ob > {output} 2> {log}"""


rule filter_offtarget_variants:
    input:
        calls=get_fixed_candidate_calls("bcf"),
        index=get_fixed_candidate_calls("bcf.csi"),
        regions="resources/target_regions/target_regions.bed",
    output:
        "results/candidate-calls/{group}.{caller}.filtered.bcf",
    params:
        extra="",
    log:
        "logs/filter_offtarget_variants/{group}.{caller}.log",
    wrapper:
        "v2.3.2/bio/bcftools/filter"


rule scanITD:
    input:
        ref=genome,
        ref_idx=genome_fai,
        regions=get_itd_regions,
        bam="results/recal/{sample}.bam",
        bai="results/recal/{sample}.bai",
    output:
        "results/candidate-calls/{sample}.ITD.vcf",
    log:
        "logs/ScanITD/{sample}.log",
    conda:
        "../envs/scanitd.yaml"
    params:
        out_name="{sample}",
        extra=config["params"]["ScanITD"]["extra"]
    threads: 2
    benchmark:
        "benchmarks/ScanITD/{sample}.tsv"
    shell:
        "(python workflow/scripts/ScanITD.py -i {input.bam} -r {input.ref} -t {input.regions} "
        "-o results/candidate-calls/{params.out_name} {params.extra}) 2> {log}"


rule make_itd_bcf:
    input:
        vcf="results/candidate-calls/{sample}.ITD.vcf",
        ref_idx=genome_fai,        
    output:
        "results/candidate-calls/{sample}.ITD.bcf"
    log:
        "logs/bcftools/reheader/ScanITD/{sample}.log",
    conda:
        "../envs/bcftools.yaml"
    params:
        itd_max_length_bp=config["params"]["ScanITD"]["itd_max_length_bp"]
    resources:
        mem_mb=8000,
    threads: 2
    shell:
        "(bcftools reheader -f {input.ref_idx} {input.vcf} | bcftools sort --max-mem {resources.mem_mb}M | "
        "bcftools view -e 'SVLEN > {params.itd_max_length_bp}' -Ob > {output}) 2> {log}"


rule bcftools_merge_ScanITD:
    input:
        calls=get_itd_bcfs,
        ind=get_itd_bcfs_index,
    output:
        "results/candidate-calls/{group}.ScanITD.bcf",
    log:
        "logs/bcf-merge/{group}.log",
    conda:
        "../envs/bcftools.yaml"
    shell:
        "(bcftools merge -m none {input.calls} | bcftools view -Ob > {output}) 2> {log}"


rule scatter_candidates:
    input:
        "results/candidate-calls/{group}.{caller}.filtered.bcf"
        if config.get("target_regions", None)
        else get_fixed_candidate_calls("bcf"),
    output:
        scatter.calling(
            "results/candidate-calls/{{group}}.{{caller}}.{scatteritem}.bcf"
        ),
    log:
        "logs/scatter-candidates/{group}.{caller}.log",
    conda:
        "../envs/rbt.yaml"
    shell:
        "rbt vcf-split {input} {output}"
