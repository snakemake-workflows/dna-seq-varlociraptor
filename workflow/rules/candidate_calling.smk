rule freebayes:
    input:
        ref="resources/genome.fasta",
        ref_idx="resources/genome.fasta.fai",
        regions=get_regions,
        # you can have a list of samples here
        samples=get_group_bams
    output:
        "results/candidate-calls/{group}.freebayes.bcf"
    log:
        "logs/freebayes/{group}.log"
    params:
        extra=config["params"].get("freebayes", ""),
    threads: 100 # use all available cores for calling
    wrapper:
        "0.55.1/bio/freebayes"


rule delly:
    input:
        ref="resources/genome.fasta",
        ref_idx="resources/genome.fasta.fai",
        samples=get_group_bams,
        index=get_group_bais,
    output:
        "results/candidate-calls/{group}.delly.bcf"
    log:
        "logs/delly/{group}.log"
    params:
        extra=config["params"].get("delly", "")
    threads: lambda _, input: len(input.samples) # delly parallelizes over the number of samples
    wrapper:
        "0.55.1/bio/delly"
