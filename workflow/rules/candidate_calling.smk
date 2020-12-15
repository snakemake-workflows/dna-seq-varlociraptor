rule freebayes:
    input:
        ref="resources/genome.fasta",
        ref_idx="resources/genome.fasta.fai",
        regions=get_regions(),
        # you can have a list of samples here
        samples=lambda w: get_group_bams(w),
        index=lambda w: get_group_bams(w, bai=True),
    output:
        "results/candidate-calls/{group}.freebayes.bcf"
    log:
        "logs/freebayes/{group}.log"
    params:
        extra="--pooled-continuous --min-alternate-count 1 {}".format(config["params"].get("freebayes", "")),
    threads: 100 # use all available cores for calling
    wrapper:
        "0.68.0/bio/freebayes"


rule delly:
    input:
        ref="resources/genome.fasta",
        ref_idx="resources/genome.fasta.fai",
        samples=lambda w: get_group_bams(w),
        index=lambda w: get_group_bams(w, bai=True),
        exclude=get_excluded_regions()
    output:
        "results/candidate-calls/{group}.delly.bcf"
    log:
        "logs/delly/{group}.log"
    params:
        extra=config["params"].get("delly", "")
    threads: lambda _, input: len(input.samples) # delly parallelizes over the number of samples
    wrapper:
        "0.68.0/bio/delly"


rule scatter_candidates:
    input:
        "results/candidate-calls/{group}.{caller}.bcf"
    output:
        scatter.calling("results/candidate-calls/{{group}}.{{caller}}.{scatteritem}.bcf")
    log:
        "logs/scatter-candidates/{group}.{caller}.log"
    conda:
        "../envs/rbt.yaml"
    shell:
        "rbt vcf-split {input} {output}"