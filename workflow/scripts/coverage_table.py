import pandas as pd

sys.stderr = open(snakemake.log[0], "w")

group_regions = dict()

for bed in [snakemake.input[0]]:
    sample = bed.split("/")[-1].split(".")[0]
    with open(bed, "r") as covered_regions:
        for line in covered_regions:
            line = line.strip().split("\t")
            chromosome = line[0]
            gene = line[3]
            coverage = line[4]
            if float(coverage) < float(snakemake.params.min_cov):
                continue
            if (chromosome, gene) not in group_regions:
                group_regions[(chromosome, gene)] = dict()
            group_regions[(chromosome, gene)][sample] = coverage

df = pd.DataFrame.from_dict(group_regions).T
df.index.names = ("chromosome", "gene")
df.reset_index(inplace=True)

with open(snakemake.output[0], "w") as csv_file:
    df.to_csv(csv_file, index=False, sep="\t")
