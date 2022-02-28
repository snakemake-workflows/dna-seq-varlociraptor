import pandas as pd


def retain_notna_cols(df):
    return df.loc[:, df.notna().sum(axis="rows")]


calls = pd.read_csv(snakemake.input[0], sep="\t")


coding = ~pd.isna(calls["hgvsp"])

retain_notna_cols(calls[~coding]).to_csv(snakemake.output.noncoding)

retain_notna_cols(calls[coding]).to_csv(snakemake.output.coding)
