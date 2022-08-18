import sys

sys.stderr = open(snakemake.log[0], "w")

import os
from pathlib import Path

import pandas as pd
import numpy as np


def join_group_hgvsgs(df):
    hgvsgs = ",".join(np.sort(pd.unique(df["hgvsg"].str.split(",").explode())))
    df.loc[:, "hgvsg"] = hgvsgs
    return df.drop_duplicates()


def join_group_consequences(df):
    consequences = ",".join(
        np.sort(pd.unique(df["consequence"].str.split(",").explode()))
    )
    df.loc[:, "consequence"] = consequences
    return df


def join_gene_vartypes(df):
    vartypes = ",".join(pd.unique(df["vartype"]))
    df = df.iloc[:1]
    df.loc[:, "vartype"] = vartypes
    return df


def load_calls(path, group):
    calls = pd.read_csv(
        path, sep="\t", usecols=["symbol", "vartype", "hgvsp", "hgvsg", "consequence"]
    )
    calls["group"] = group
    calls.loc[:, "consequence"] = calls["consequence"].str.replace("&", ",")
    return calls.drop_duplicates()


def load_group_annotation():
    if "group_annotation" in snakemake.input.keys():
        group_annotation = (
            pd.read_csv(
                snakemake.input.group_annotation, sep="\t", dtype={"group": str}
            )
            .set_index("group")
            .sort_index()
        )
        group_annotation.loc[
            set(snakemake.params.groups) - set(group_annotation.index)
        ] = [pd.NA] * len(group_annotation.columns)
    else:
        group_annotation = pd.DataFrame({"group": snakemake.params.groups}).set_index(
            "group"
        )
    group_annotation = group_annotation.T
    group_annotation.index.name = "annotation"
    return group_annotation


def sort_by_recurrence(matrix, no_occurence_check_func):
    matrix["nocalls"] = no_occurence_check_func(matrix).sum(axis="columns")
    matrix = matrix.sort_values("nocalls", ascending=True).drop(
        labels=["nocalls"], axis="columns"
    )
    return matrix


def add_missing_groups(matrix, groups, index_mate):
    for group in groups:
        if (index_mate, group) not in matrix.columns:
            matrix[(index_mate, group)] = pd.NA
    return matrix


def attach_group_annotation(matrix, group_annotation):
    index_cols = matrix.index.names
    return pd.concat(
        [group_annotation.reset_index(drop=True), matrix.reset_index()]
    ).set_index(index_cols)


def gene_oncoprint(calls, group_annotation):
    calls = calls[["group", "symbol", "vartype", "consequence"]]
    if not calls.empty:
        grouped = (
            calls.drop_duplicates().groupby(["symbol"]).apply(join_group_consequences)
        )
        grouped = (
            grouped.drop_duplicates()
            .groupby(["group", "symbol"])
            .apply(join_gene_vartypes)
        )
        matrix = grouped.set_index(["symbol", "consequence", "group"]).unstack(
            level="group"
        )
        matrix = add_missing_groups(matrix, snakemake.params.groups, "vartype")
        matrix.columns = matrix.columns.droplevel(0)  # remove superfluous header
        if len(matrix.columns) > 1:
            # sort by recurrence
            matrix = sort_by_recurrence(matrix, lambda matrix: matrix.isna())

        matrix = attach_group_annotation(matrix, group_annotation)
        return matrix.reset_index()
    else:
        cols = ["symbol", "consequence"] + snakemake.params.groups
        return pd.DataFrame({col: [] for col in cols})


def variant_oncoprint(gene_calls, group_annotation):
    gene_calls = gene_calls[["group", "hgvsp", "hgvsg", "consequence"]]
    gene_calls.loc[:, "exists"] = "X"
    grouped = gene_calls.drop_duplicates().groupby(["hgvsp"]).apply(join_group_hgvsgs)
    matrix = grouped.set_index(["hgvsp", "hgvsg", "consequence", "group"]).unstack(
        level="group"
    )

    matrix = add_missing_groups(matrix, snakemake.params.groups, "exists")
    matrix.columns = matrix.columns.droplevel(0)  # remove superfluous header

    if len(matrix.columns) > 1:
        # sort by recurrence
        matrix = sort_by_recurrence(matrix, lambda matrix: matrix.isna())

    matrix = attach_group_annotation(matrix, group_annotation)

    return matrix.reset_index()


calls = pd.concat(
    [
        load_calls(path, sample)
        for path, sample in zip(snakemake.input.calls, snakemake.params.groups)
    ]
)

group_annotation = load_group_annotation()

gene_oncoprint(calls, group_annotation).to_csv(
    snakemake.output.gene_oncoprint, sep="\t", index=False
)

os.makedirs(snakemake.output.variant_oncoprints)
for gene, gene_calls in calls.groupby("symbol"):
    variant_oncoprint(gene_calls, group_annotation).to_csv(
        Path(snakemake.output.variant_oncoprints) / f"{gene}.tsv", sep="\t", index=False
    )
