import sys

sys.stderr = open(snakemake.log[0], "w")

import os
from pathlib import Path

import pandas as pd
import numpy as np


def join_gene_variants(df):
    vartypes = ",".join(pd.unique(df["vartype"]))
    consequences = ",".join(np.sort(pd.unique(df["consequence"].str.split(",").explode())))
    df = df.iloc[:1]
    df.loc[:, "vartype"] = vartypes
    df.loc[:, "consequence"] = consequences
    return df


def load_calls(path, group):
    calls = pd.read_csv(path, sep="\t", usecols=["symbol", "vartype", "hgvsp", "hgvsg", "consequence"])
    calls["group"] = group
    calls.loc[:, "consequence"] = calls["consequence"].str.replace("&", ",")
    return calls.drop_duplicates()


def gene_oncoprint(calls):
    matrix = (
        calls[["group", "symbol", "vartype", "consequence"]].drop_duplicates().groupby(["group", "symbol"]).apply(join_gene_variants)
        .set_index(["symbol", "consequence", "group"])
        .unstack(level="group")
    )
    matrix.columns = matrix.columns.droplevel(0)  # remove superfluous header
    matrix["nocalls"] = matrix.isna().sum(axis="columns")
    return matrix.sort_values("nocalls", ascending=False).drop(labels=["nocalls"], axis="columns").reset_index()


def variant_oncoprint(gene_calls):
    gene_calls = gene_calls[["group", "hgvsp", "hgvsg", "consequence"]]
    gene_calls.loc[:, "exists"] = 1
    matrix = (
        gene_calls.set_index(["hgvsp", "hgvsg", "consequence", "group"])
        .unstack(level="group")
        .fillna(0)
    )
    matrix.columns = matrix.columns.droplevel(0)  # remove superfluous header

    matrix["nocalls"] = (matrix == 0).sum(axis="columns")
    return matrix.sort_values("nocalls", ascending=False).drop(labels=["nocalls"], axis="columns").reset_index()


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
