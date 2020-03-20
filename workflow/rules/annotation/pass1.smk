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
        "71dc35eb2b0287ae3b3633b35c9af49b780d1850/bio/snpeff"
