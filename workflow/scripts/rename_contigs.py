import pysam

sys.stderr = open(snakemake.log[0], "w")

rename_expressions = snakemake.params.expressions

vcf_file = pysam.VariantFile(snakemake.input[0])
contigs = list(vcf_file.header.contigs)
vcf_file.close()

with open(snakemake.output[0], "w") as out:
    for contig in contigs:
        for expression in snakemake.params.expressions:
            updated_contig = eval(expression)
        print(f"{contig}\t{updated_contig}", file=out)
