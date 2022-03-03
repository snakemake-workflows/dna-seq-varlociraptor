import pandas as pd

sys.stderr = open(snakemake.log[0], "w")

calls = pd.read_csv(snakemake.input[0], sep="\t")


coding = ~pd.isna(calls["hgvsp"])

calls[~coding].dropna(how="all").to_csv(snakemake.output.noncoding, index=False)

calls[coding].dropna(how="all").to_csv(snakemake.output.coding, index=False)