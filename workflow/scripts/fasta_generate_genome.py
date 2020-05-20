nbr_chromosomes = snakemake.input[1]
out_file = open(snakemake.output[0], "w")
with open(snakemake.input[0], "r") as ref_idx:
    for i in range(nbr_chromosomes):
        line = ref_idx.readline()
        fields = line.strip().split("\t")
        chromosome = fields[0]
        chromosome_len = fields[1]
        print("chr{chr}\t{len}".format(chr=chromosome, len=chromosome_len), file=out_file)
out_file.close()
