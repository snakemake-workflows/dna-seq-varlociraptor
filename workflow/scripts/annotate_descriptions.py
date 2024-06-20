import pandas as pd

sys.stderr = open(snakemake.log[0], "w")

signatures_df = pd.read_csv(snakemake.input.sig, sep="\t")
description_df = pd.read_csv(snakemake.input.desc, sep="\t")
signatures_df = pd.merge(signatures_df, description_df, how="left", on="Signature")

signatures_df.to_csv(snakemake.output[0], sep="\t", index=False)
