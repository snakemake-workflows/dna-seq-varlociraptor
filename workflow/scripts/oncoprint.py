import sys

sys.stderr = open(snakemake.log[0], "w")

import os
from pathlib import Path

import pandas as pd
import numpy as np


def join_gene_variants(df):
    vartypes = ",".join(pd.unique(df["vartype"]))
    consequences = ",".join(
        np.sort(pd.unique(df["consequence"].str.split(",").explode()))
    )
    df = df.iloc[:1]
    df.loc[:, "vartype"] = vartypes
    df.loc[:, "consequence"] = consequences
    return df


def load_calls(path, group):
    calls = pd.read_csv(
        path, sep="\t", usecols=["symbol", "vartype", "hgvsp", "hgvsg", "consequence"]
    )
    calls["group"] = group
    calls.loc[:, "consequence"] = calls["consequence"].str.replace("&", ",")
    return calls.drop_duplicates()


def sort_by_recurrence(matrix, no_occurence_check_func):
    matrix["nocalls"] = no_occurence_check_func(matrix).sum(axis="columns")
    matrix = matrix.sort_values("nocalls", ascending=False).drop(
        labels=["nocalls"], axis="columns"
    )
    return matrix


def add_missing_groups(matrix, groups, index_mate):
    for group in groups:
        if (index_mate, group) not in matrix.columns:
            matrix[(index_mate, group)] = 0
    return matrix


def gene_oncoprint(calls):
    calls = calls[["group", "symbol", "vartype", "consequence"]]
    if not calls.empty:
        grouped = (
            calls.drop_duplicates()
            .groupby(["group", "symbol"])
            .apply(join_gene_variants)
        )
        matrix = grouped.set_index(["symbol", "consequence", "group"]).unstack(
            level="group"
        )
        matrix = add_missing_groups(matrix, snakemake.params.groups, "vartype")
        matrix.columns = matrix.columns.droplevel(0)  # remove superfluous header
        if len(matrix.columns) > 1:
            # sort by recurrence
            matrix = sort_by_recurrence(matrix, lambda matrix: matrix.isna())
        return matrix.reset_index()
    else:
        return calls


def variant_oncoprint(gene_calls):
    gene_calls = gene_calls[["group", "hgvsp", "hgvsg", "consequence"]]
    gene_calls.loc[:, "exists"] = 1
    matrix = (
        gene_calls.set_index(["hgvsp", "hgvsg", "consequence", "group"])
        .unstack(level="group")
        .fillna(0)
    )

    matrix = add_missing_groups(matrix, snakemake.params.groups, "exists")
    matrix.columns = matrix.columns.droplevel(0)  # remove superfluous header

    if len(matrix.columns) > 1:
        # sort by recurrence
        matrix = sort_by_recurrence(matrix, lambda matrix: matrix == 0)
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
