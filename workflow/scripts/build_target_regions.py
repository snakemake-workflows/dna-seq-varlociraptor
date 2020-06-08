log_file=open(snakemake.log[0], 'w')
out=open(snakemake.output[0], 'w')

with open(snakemake.input[0], 'r') as primer_file:
    for line in primer_file.readlines():
        line = line.strip().split("\t")
        chr_1 = line[0]
        chr_2 = line[1]
        if chr_1 == chr_2:
            print("{chr}\t{start}\t{end}".format(chr=chr_1, start=line[1], end=line[5]), file=out)
        else:
            print("{primer}\t{chr1}\t{chr2}".format(primer=line[6],chr1=chr_1, chr2=chr_2), file=log_file)