rule download_snpeff_db:
    output:
        directory("{working_dir}/results/databases/{{ref}}".format(working_dir=os.getcwd()))
    params:
        db_dir=os.path.join(os.getcwd(), "results/databases/")
    cache: True
    conda:
        "../../envs/snpeff.yaml"
    shell:
        "snpEff download -dataDir {params.db_dir} {wildcards.ref}"

rule snpeff:
    input:
        calls="results/calls/{group}.bcf",
        db="{working_dir}/results/databases/{build}.{release}".format(working_dir=os.getcwd(), build=config["ref"]["build"], release=config["ref"]["snpeff_release"])
    output:
        calls="results/calls/{group}.annotated.bcf",
        stats="results/snpeff/{group}.html",
        csvstats="results/snpeff/{group}.csv"
    log:
        "logs/snpeff/{group}.log"
    params:
        reference="{}.{}".format(config["ref"]["build"], config["ref"]["snpeff_release"]),
        data_dir = os.path.join(os.getcwd(), "results/databases/"),
        extra="-Xmx4g -nodownload"
    wrapper:
        "0.50.4/bio/snpeff"
