import sys
import pandas as pd

sys.stderr = open(snakemake.log[0], "w")


def write(df, path):
    df.drop(["canonical"], axis="columns").dropna(how="all", axis="columns").to_csv(
        path, index=False, sep=","
    )


calls = pd.read_csv(snakemake.input[0], sep="\t")
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
