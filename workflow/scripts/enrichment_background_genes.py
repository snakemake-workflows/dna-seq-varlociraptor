import sys
import pandas as pd

sys.stderr = open(snakemake.log[0], "w")

df = pd.read_csv(snakemake.input.coverage, sep="\t")

sample_cols = [c for c in df.columns if c not in ("chromosome", "gene")]
df[sample_cols] = df[sample_cols].apply(pd.to_numeric, errors="coerce")
df["mean_coverage"] = df[sample_cols].mean(axis=1)

background = df.loc[df["mean_coverage"] >= snakemake.params.min_cov, "gene"]

with open(snakemake.output[0], "w") as f:
    for gene in background:
        f.write(f"{gene}\n")

print(f"Background genes written: {len(background)}", file=sys.stderr)
