import sys

sys.stderr = open(snakemake.log[0], "w")

import os
from pathlib import Path
import json

import pandas as pd
import numpy as np

from sklearn.feature_selection import chi2
from statsmodels.stats.multitest import fdrcorrection


def join_group_hgvsgs(df):
    hgvsgs = ",".join(np.sort(pd.unique(df["HGVSg"].str.split(",").explode())))
    df.loc[:, "HGVSg"] = hgvsgs
    return df.drop_duplicates()


def join_group_consequences(df):
    consequences = ",".join(
        np.sort(pd.unique(df["Consequence"].str.split(",").explode()))
    )
    df.loc[:, "Consequence"] = consequences
    return df


def join_gene_vartypes(df):
    vartypes = ",".join(pd.unique(df["vartype"]))
    df = df.iloc[:1]
    df.loc[:, "vartype"] = vartypes
    return df


def load_calls(path, group):
    calls = pd.read_csv(path, sep="\t")
    calls = calls[["Symbol", "vartype", "HGVSp", "HGVSg", "Consequence"] + get_allelefreq_columns(calls)]
    calls["group"] = group
    calls.loc[:, "Consequence"] = calls["Consequence"].str.replace("&", ",")
    return calls.drop_duplicates()


def get_allelefreq_columns(calls):
    return list(calls.columns[calls.columns.str.endswith("allele frequency")].values)


def load_group_annotation():
    if snakemake.input.group_annotation:
        group_annotation = (
            pd.read_csv(
                snakemake.input.group_annotation, sep="\t", dtype={"group": str}
            )
            .set_index("group")
            .sort_index()
        )
        group_annotation = group_annotation.loc[snakemake.params.groups]
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
    return (
        pd.concat([group_annotation.reset_index(drop=True), matrix.reset_index()])
        .set_index(index_cols)
        .reset_index()
    )


def gene_by_sample_oncoprint(calls):
    calls = calls[["group", "Symbol", "vartype", "Consequence"]]
    if not calls.empty:
        grouped = (
            calls.drop_duplicates().groupby(["Symbol"]).apply(join_group_consequences)
        ).reset_index(drop=True)
        grouped = (
            grouped.drop_duplicates()
            .groupby(["group", "Symbol"])
            .apply(join_gene_vartypes)
        )
        matrix = grouped.set_index(["Symbol", "Consequence", "group"]).unstack(
            level="group"
        )
        matrix = add_missing_groups(matrix, snakemake.params.groups, "vartype")
        matrix.columns = matrix.columns.droplevel(0)  # remove superfluous header
        if len(matrix.columns) > 1:
            # sort by recurrence
            matrix = sort_by_recurrence(matrix, lambda matrix: matrix.isna())

        return matrix
    else:
        cols = ["Symbol", "Consequence"] + list(snakemake.params.groups)
        return pd.DataFrame({col: [] for col in cols}).set_index(
            list(snakemake.params.groups)
        )


def group_by_variant_oncoprint(calls, group_annotation):
    allelefreq_cols = get_allelefreq_columns(calls)
    calls = calls[["group"] + snakemake.params.group_by_variant_oncoprint_header_labels + allelefreq_cols]

    def join_allele_freqs(row):
        return " + ".join(colname.split(":")[0].strip() for colname, vaf in row[allelefreq_cols].iteritems() if vaf > 0.0)

    calls["status"] = calls.apply(join_allele_freqs, axis="columns")
    calls.drop(columns=allelefreq_cols, inplace=True)
    # unstack and drop the status column (with label level_0 after reset_index)
    calls = calls.set_index(["group"] + snakemake.params.group_by_variant_oncoprint_header_labels).unstack(level="group").T.reset_index().drop(columns=["level_0"]).set_index("group")
    # backup column index as it will be destroyed by the pd.concat below
    colindex = calls.columns
    # add group annotation and set as index
    calls = pd.concat([group_annotation.T, calls]).set_index(group_annotation.columns.tolist(), append=True)
    # restore column index
    calls.columns = colindex
    # reset index so that table has the right format (no column labels)
    calls.reset_index()
    return calls



def variant_oncoprint(gene_calls, group_annotation):
    gene_calls = gene_calls[["group", "HGVSp", "HGVSg", "Consequence"]]
    gene_calls.loc[:, "exists"] = "X"
    grouped = gene_calls.drop_duplicates().groupby(["HGVSp"]).apply(join_group_hgvsgs)
    matrix = grouped.set_index(["HGVSp", "HGVSg", "Consequence", "group"]).unstack(
        level="group"
    )

    matrix = add_missing_groups(matrix, snakemake.params.groups, "exists")
    matrix.columns = matrix.columns.droplevel(0)  # remove superfluous header

    if len(matrix.columns) > 1:
        # sort by recurrence
        matrix = sort_by_recurrence(matrix, lambda matrix: matrix.isna())

    matrix = attach_group_annotation(matrix, group_annotation)

    return matrix


def store(data, output, labels_df, label_idx=None):
    _labels_df = labels_df
    if label_idx is not None:
        _labels_df = labels_df.iloc[[label_idx]]

    # add labels
    index_cols = data.index.names
    cols = data.columns
    data = pd.concat([_labels_df, data.reset_index()]).set_index(index_cols)
    # restore column order
    data = data[cols]

    data.to_csv(output, sep="\t")


def sort_oncoprint_labels(data):
    labels_df = snakemake.params.labels
    labels = labels_df.index

    for label_idx, label in enumerate(labels):
        outdata = data
        if not data.empty:
            feature_matrix = data.reset_index(drop=True).T.copy()
            feature_matrix[~pd.isna(feature_matrix)] = True
            feature_matrix[pd.isna(feature_matrix)] = False

            # target vector: label values, converted into factors
            target_vector = labels_df.loc[label]
            # ignore any NA in the target vector and correspondingly remove the rows in the feature matrix
            not_na_target_vector = target_vector[~pd.isna(target_vector)].astype("category")
            feature_matrix = feature_matrix.loc[not_na_target_vector.index]

            # calculate mutual information for 100 times and take the mean for each feature
            _, pvals = chi2(feature_matrix, not_na_target_vector)
            sorted_idx = np.argsort(pvals)

            _, fdr = fdrcorrection(pvals)

            # clone data
            sorted_data = data.copy(deep=True)

            # sort by label
            sorted_target_vector = target_vector.sort_values()
            sorted_data = sorted_data[sorted_target_vector.index]

            # add mutual information
            sorted_data.insert(0, "FDR dependency", np.around(fdr, 3))
            sorted_data.insert(0, "p-value dependency", np.around(pvals, 3))

            outdata = sorted_data.iloc[sorted_idx]
        outpath = os.path.join(snakemake.output.gene_oncoprint_sortings, f"{label}.tsv")
        store(outdata, outpath, labels_df, label_idx=label_idx)


calls = pd.concat(
    [
        load_calls(path, sample)
        for path, sample in zip(snakemake.input.calls, snakemake.params.groups)
    ]
)

gene_oncoprint = gene_by_sample_oncoprint(calls)

group_annotation = load_group_annotation()
gene_oncoprint_main = attach_group_annotation(gene_oncoprint, group_annotation)
gene_oncoprint_main.to_csv(snakemake.output.gene_oncoprint, sep="\t", index=False)

group_by_variant_oncoprint = group_by_variant_oncoprint(calls, group_annotation)
group_by_variant_oncoprint.to_csv(snakemake.output.group_by_variant_oncoprint, sep="\t", index=False)

os.makedirs(snakemake.output.gene_oncoprint_sortings)

sort_oncoprint_labels(gene_oncoprint)


os.makedirs(snakemake.output.variant_oncoprints)
for gene, gene_calls in calls.groupby("Symbol"):
    variant_oncoprint(gene_calls, group_annotation).to_csv(
        Path(snakemake.output.variant_oncoprints) / f"{gene}.tsv", sep="\t", index=False
    )
