import re
import sys
import numpy as np
import pandas as pd

sys.stderr = open(snakemake.log[0], "w")


def write(df, path):
    df.drop(["canonical"], axis="columns").dropna(how="all", axis="columns").to_csv(
        path, index=False, sep="\t"
    )


def format_floats(df):
    for col_name in df:
        if issubclass(df[col_name].dtype.type, np.floating):
            df[col_name] = [
                "{:.2e}".format(x) if x < 0.1 and x > 0 else round(x, 2)
                for x in df[col_name]
            ]
    return df


def trim_hgvsp_entries(df):
    if "hgvsp" in df.columns:
        df["hgvsp"] = df["hgvsp"].apply(lambda x: re.sub(r"^.*?:", "", x) if pd.notnull(x) else x)
    return df


calls = pd.read_csv(snakemake.input[0], sep="\t")
calls = format_floats(calls)
calls = trim_hgvsp_entries(calls)
calls.set_index("gene", inplace=True, drop=False)


coding = ~pd.isna(calls["hgvsp"])
canonical = calls["canonical"]

write(calls[~coding & canonical], snakemake.output.noncoding)

# coding variants
# Here we have all variant info in hgvsp,
# hence we can drop the other variant specific columns.
write(
    calls[coding & canonical].drop(
        [
            "chromsome",
            "position",
            "reference allele",
            "alternative allele",
            "end position",
            "event",
            "id",
        ],
        axis="columns",
    ),
    snakemake.output.coding,
)
# TODO add possibility to also see non-canonical transcripts (low priority, once everything else works).
