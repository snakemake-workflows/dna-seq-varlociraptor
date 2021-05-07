import sys
import pandas as pd

sys.stderr = open(snakemake.log[0], "w")

data = pd.read_csv(snakemake.input.tsv, sep="\t")
data.to_excel(snakemake.output.xlsx, index=False)
