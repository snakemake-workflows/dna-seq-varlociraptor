import sys
import pandas as pd

sys.stderr = open(snakemake.log[0], "w")

genes = set()
for tsv in snakemake.input:
    df = pd.read_csv(tsv, sep="\t")
    if "SYMBOL" in df.columns:
        genes.update(df["SYMBOL"].dropna().unique())

with open(snakemake.output[0], "w") as f:
    for gene in sorted(genes):
        f.write(f"{gene}\n")

print(f"Gene list written: {len(genes)} genes", file=sys.stderr)
