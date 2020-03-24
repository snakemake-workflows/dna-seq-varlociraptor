rule freebayes:
    threads:
        100 # use all available cores for calling
    input:
        ref="results/refs/genome.fasta",
        ref_idx="results/refs/genome.fasta.fai",
        # you can have a list of samples here
        samples=get_group_bams
    output:
        "results/candidate-calls/{group}.freebayes.bcf"
    log:
        "logs/freebayes/{group}.log"
    params:
        extra=config["params"].get("freebayes", ""),
        chunksize=100000
    wrapper:
        "0.50.3/bio/freebayes"

rule delly:
    threads:
        4 # this should be the number of samples in the group, since delly is parallized over the samples
    input:
        ref="results/refs/genome.fasta",
        ref_idx="results/refs/genome.fasta.fai",
        samples=get_group_bams,
        index=get_group_bais,
    output:
        "results/candidate-calls/{group}.delly.bcf"
    params:
        extra=config["params"].get("delly", "")
    log:
        "logs/delly/{group}.log"
    wrapper:
        "0.43.0/bio/delly"
