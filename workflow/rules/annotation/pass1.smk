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
        "2421904a0c01588f5c67d3e7f68024ccc9f6787b/bio/snpeff"
