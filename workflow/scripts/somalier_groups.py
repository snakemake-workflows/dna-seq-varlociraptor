import sys
sys.stderr = open(snakemake.log[0], "w")

snakemake.params.samples.groupby("group").apply(
    lambda df: ",".join(df["sample_name"]), include_groups=False
).to_csv(snakemake.output[0], index=False, header=False, sep="\t")