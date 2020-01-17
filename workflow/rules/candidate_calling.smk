rule freebayes:
    input:
        ref="refs/genome.fasta",
        ref_idx="refs/genome.fasta.fai",
        # you can have a list of samples here
        samples=get_group_bams
    output:
        "candidate-calls/{group}.freebayes.bcf"
    log:
        "logs/freebayes/{group}.log"
    params:
        extra=config["params"].get("freebayes", ""),
        chunksize=100000
    wrapper:
        "0.43.0/bio/freebayes"

rule delly:
    input:
        ref="refs/genome.fasta",
        ref_idx="refs/genome.fasta.fai",
        samples=get_group_bams,
        index=get_group_bais,
    output:
        "candidate-calls/{group}.delly.bcf"
    params:
        extra=config["params"].get("delly", "")
    log:
        "logs/delly/{group}.log"
    wrapper:
        "0.43.0/bio/delly"
