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
    df.to_csv(path, index=False, sep="\t")


calls = pd.read_csv(snakemake.input[0], sep="\t", usecols=["chromosome", "position", "id", "mateid", "feature_id", "feature_name", "exon"], dtype={"exon": "Int64"})
calls[["feature_name", "feature_id"]] = calls[["feature_id", "feature_name"]].fillna("('intronic',)")

# Obtain first entry of columns annotated as tuple by arriba
for col in ["feature_id", "feature_name"]:
    calls[col] = (
        calls[col]
        .apply(eval)
        .apply(lambda x: x[0])
    )

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