rule snpeff:
    input:
        "results/calls/{group}.bcf"
    output:
        calls="results/calls/{group}.annotated.bcf",
        stats="results/snpeff/{group}.html",
        csvstats="results/snpeff/{group}.csv"
    log:
        "logs/snpeff/{group}.log"
    params:
        reference="{}.{}".format(config["ref"]["build"], config["ref"]["snpeff_release"]),
        extra="-Xmx4g"
    wrapper:
        "0.50.4/bio/snpeff"
