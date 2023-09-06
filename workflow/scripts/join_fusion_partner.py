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


calls = pd.read_csv(snakemake.input[0], sep="\t")

calls = calls[["chromosome", "position", "id", "mateid", "symbol", "exon"]]
calls = calls.dropna(how="any", axis="rows")
calls = calls.drop_duplicates()
first_fusions = calls["id"].str.contains("a$")
first_calls = calls[first_fusions]
second_calls = calls[~first_fusions]
paired_fusions = first_calls.merge(
    second_calls,
    left_on="mateid",
    right_on="id",
    suffixes=(" partner1", " partner2"),
)
paired_fusions = paired_fusions.filter(regex="^(?!mateid|id)")

write(paired_fusions, snakemake.output.fusions)