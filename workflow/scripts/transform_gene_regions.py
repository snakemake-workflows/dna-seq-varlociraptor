sys.stderr = open(snakemake.log[0], "w")

out_file = open(snakemake.output[0], "w")

with open(snakemake.input[0], "r") as annotations:
    for line in annotations.readlines():
        if line.startswith("#"):
            continue
        line = line.split("\t")
        chromosome = line[0]
        feature = line[2]
        if feature != "gene" or chromosome not in list(map(str, range(1, 23))) + [
            "X",
            "Y",
        ]:
            continue
        start = str(int(line[3]) - 1)
        end = line[4]
        desc = dict([x.split(" ") for x in line[8].split("; ")])
        gene = desc.get("gene_name", desc["gene_id"])[1:-1]
        print("\t".join([chromosome, start, end, gene]), file=out_file)
out_file.close()
