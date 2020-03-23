rule snpeff:
    input:
        "calls/{group}.bcf"
    output:
        calls="calls/{group}.annotated.bcf",
        stats="snpeff/{group}.html",
        csvstats="snpeff/{group}.csv"
    log:
        "logs/snpeff/{group}.log"
    params:
        reference="{}.{}".format(config["ref"]["build"], config["ref"]["snpeff_release"]),
        extra="-Xmx4g"
    wrapper:
        "0.50.4/bio/snpeff"
