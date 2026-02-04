rule freebayes:
    input:
        ref=access.random(genome),
        ref_idx=genome_fai,
        regions="results/regions/{group}.expanded_regions.filtered.bed",
        # you can have a list of samples here
        alns=access.random(lambda w: get_group_bams(w)),
        idxs=lambda w: get_group_bams(w, crai=True),
    output:
        "results/candidate-calls/freebayes/{group}/{group}.bcf",
    log:
        "logs/freebayes/{group}.log",
    params:
        # genotyping is performed by varlociraptor, hence we deactivate it in freebayes by
        # always setting --pooled-continuous
        extra="--pooled-continuous --min-alternate-count {} --min-alternate-fraction {} {}".format(
            1 if is_activated("calc_consensus_reads") else 2,
            config["params"]["freebayes"].get("min_alternate_fraction", "0.05"),
            config["params"]["freebayes"].get("extra", ""),
        ),
    threads: max(workflow.cores - 1, 1)  # use all available cores -1 (because of the pipe) for calling
    wrapper:
        "v8.0.3/bio/freebayes"


rule delly:
    input:
        ref=access.random(genome),
        ref_idx=genome_fai,
        alns=access.random(lambda w: get_group_bams(w)),
        index=lambda w: get_group_bams(w, crai=True),
        exclude=get_delly_excluded_regions(),
    output:
        "results/candidate-calls/delly/{group}/{group}.bcf",
    log:
        "logs/delly/{group}.log",
    params:
        extra=config["params"].get("delly", ""),
    threads: lambda _, input: len(input.alns)  # delly parallelizes over the number of samples
    wrapper:
        "v8.0.3/bio/delly"


# Delly breakends lead to invalid BCFs after VEP annotation (invalid RLEN). Therefore we exclude them for now.
rule fix_delly_calls:
    input:
        "results/candidate-calls/delly/{group}/{group}.bcf",
    output:
        "results/candidate-calls/delly/{group}/{group}.no_bnds.bcf",
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
        "results/candidate-calls/{caller}/filtered/{group}/{group}.bcf",
    params:
        extra="",
    log:
        "logs/filter_offtarget_variants/{group}/{group}.{caller}.log",
    wrapper:
        "v2.3.2/bio/bcftools/filter"


rule scatter_candidates:
    input:
        "results/candidate-calls/{caller}/filtered/{group}/{group}.bcf"
        if config.get("target_regions", None)
        else get_fixed_candidate_calls("bcf"),
    output:
        scatter.calling(
            "results/candidate-calls/{{caller}}/{{group}}/{{group}}.{scatteritem}.bcf"
        ),
    log:
        "logs/scatter-candidates/{caller}/{group}/{group}.log",
    conda:
        "../envs/rbt.yaml"
    shell:
        "rbt vcf-split {input} {output}"
