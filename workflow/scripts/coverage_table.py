import pandas as pd
import numpy as np

sys.stderr = open(snakemake.log[0], "w")


def add_missing_columns(df, samples):
    missing_columns = set(samples).difference(df.columns)
    df[list(missing_columns)] = np.nan
    return df


group_regions = dict()
samples = []
for bed in snakemake.input:
    sample = bed.split("/")[-1].split(".")[0]
    samples.append(sample)
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

if bool(group_regions):
    df = pd.DataFrame.from_dict(group_regions).T
    df.index.names = ("chromosome", "gene")
    df.reset_index(inplace=True)
    df = add_missing_columns(df, samples)
else:
    df = pd.DataFrame(columns=["chromosome", "gene"] + samples)

with open(snakemake.output[0], "w") as csv_file:
    df.to_csv(csv_file, index=False, sep="\t")
