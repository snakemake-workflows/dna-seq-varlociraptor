import sys

sys.stderr = open(snakemake.log[0], "w")

import os
from pathlib import Path

import pandas as pd


def join_vartypes(df):
    df["vartype"] = ",".join(df["vartype"])
    return df.drop_duplicates()


def load_calls(path, group):
    calls = pd.read_csv(path, sep="\t", usecols=["symbol", "vartype", "hgvsp", "hgvsg"])
    calls["group"] = group
    return calls.drop_duplicates()


def gene_oncoprint(calls):
    matrix = (
        calls[["group", "symbol", "vartype"]]
        .drop_duplicates()
        .groupby(["group", "symbol"])
        .apply(join_vartypes)
        .set_index(["symbol", "group"])
        .unstack(level="group")
    )
    matrix.columns = matrix.columns.droplevel(0)  # remove superfluous header
    return matrix.reset_index()


def variant_oncoprint(gene_calls):
    gene_calls = gene_calls[["group", "hgvsp", "hgvsg"]]
    gene_calls.loc[:, "exists"] = 1
    matrix = (
        gene_calls.set_index(["hgvsp", "hgvsg", "group"])
        .unstack(level="group")
        .fillna(0)
    )
    matrix.columns = matrix.columns.droplevel(0)  # remove superfluous header

    return matrix.reset_index()


calls = pd.concat(
    [
        load_calls(path, sample)
        for path, sample in zip(snakemake.input, snakemake.params.groups)
    ]
)


gene_oncoprint(calls).to_csv(snakemake.output.gene_oncoprint, sep="\t", index=False)

os.makedirs(snakemake.output.variant_oncoprints)
for gene, gene_calls in calls.groupby("symbol"):
    variant_oncoprint(gene_calls).to_csv(
        Path(snakemake.output.variant_oncoprints) / f"{gene}.tsv", sep="\t", index=False
    )
