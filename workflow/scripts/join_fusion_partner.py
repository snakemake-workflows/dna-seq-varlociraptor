import sys

sys.stderr = open(snakemake.log[0], "w")

import json
import re

import numpy as np
import pandas as pd

def write(df, path):
    if not df.empty:
        remaining_columns = df.dropna(how="all", axis="columns").columns.tolist()
        df = df[remaining_columns]
        canonical_columns = df.filter(regex="^canonical_partner", axis=1).columns.tolist()
        df = df.drop(labels=canonical_columns, axis=1)
    df.to_csv(path, index=False, sep="\t")


calls = pd.read_csv(snakemake.input[0], sep="\t", dtype={"exon": "Int64"})

calls = calls[["chromosome", "position", "id", "mateid", "symbol", "exon", "canonical"]]
calls = calls.drop_duplicates()

canonical = calls["canonical"].notnull()
calls = calls[canonical]

first_partners = calls["id"].str.contains("a$")
first_calls = calls[first_partners]
second_calls = calls[~first_partners]
paired_fusions = first_calls.merge(
    second_calls,
    left_on="mateid",
    right_on="id",
    suffixes=("_partner1", "_partner2"),
)
paired_fusions = paired_fusions.filter(regex="^(?!mateid|id)")

write(paired_fusions, snakemake.output.fusions)