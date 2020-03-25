rule download_snpeff_db:
    output:
        directory("results/refs/snpeff/{ref}")
    params:
        db_dir=lambda _, output: Path(output[0]).parent.resolve()
    cache: True
    conda:
        "../../envs/snpeff.yaml"
    shell:
        "snpEff download -dataDir {params.db_dir} {wildcards.ref}"

rule snpeff:
    input:
        calls="results/calls/{group}.bcf",
        db="results/refs/snpeff/{build}.{snpeff_release}".format(**config["ref"])
    output:
        calls="results/calls/{group}.annotated.bcf",
        stats="results/snpeff/{group}.html",
        csvstats="results/snpeff/{group}.csv"
    log:
        "logs/snpeff/{group}.log"
    params:
        reference="{build}.{snpeff_release}".format(**config["ref"]),
        data_dir=lambda _, input: Path(input.db).parent.resolve(),
        extra="-Xmx4g -nodownload"
    resources:
        mem_mb=4000
    wrapper:
        "0.50.4/bio/snpeff"
