import pandas as pd


calls = pd.read_csv(snakemake.input[0], sep="\t")


coding = ~pd.isna(calls["hgvsp"])

calls[~coding].dropna(how="all").to_csv(snakemake.output.noncoding, index=False)

calls[coding].dropna(how="all").to_csv(snakemake.output.coding, index=False)
