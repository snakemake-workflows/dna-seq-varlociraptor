import sys

sys.stderr = open(snakemake.log[0], "w")

import json
import re

import numpy as np
import pandas as pd


PROB_EPSILON = 0.01  # columns with all probabilities below will be dropped


def write(df, path):
    df = df.drop(["canonical"], axis="columns", errors="ignore")
    if not df.empty:
        remaining_columns = df.dropna(how="all", axis="columns").columns.tolist()
        if path == snakemake.output.coding:
            # ensure that these columns are kept, even if they contain only NAs in a coding setting
            remaining_columns.extend(["revel", "hgvsp", "symbol"])
            remaining_columns = [ col for col in df.columns if col in remaining_columns ]
        df = df[remaining_columns]
    df.to_csv(path, index=False, sep="\t")


def format_floats(df):
    for col_name in df:
        if issubclass(df[col_name].dtype.type, np.floating):
            df[col_name] = [
                "{:.2e}".format(x) if x < 0.1 and x > 0 else round(x, 2)
                for x in df[col_name]
            ]
    return df


def drop_cols_by_predicate(df, columns, predicate):
    predicate_true_cols = [
        col for col in columns if predicate(df[col].astype(float)).all()
    ]
    return df.drop(labels=predicate_true_cols, axis="columns")


def drop_low_prob_cols(df):
    return drop_cols_by_predicate(
        df,
        df.columns[df.columns.str.startswith("prob: ")],
        lambda probs: probs <= PROB_EPSILON,
    )


def get_vaf_columns(df):
    return df.columns[df.columns.str.endswith(": allele frequency")]


def get_samples(df):
    return list(get_vaf_columns(df).str.replace(": allele frequency", ""))


def get_vartype(rec):
    ref_allele = rec["reference allele"]
    alt_allele = rec["alternative allele"]
    if alt_allele == "<DEL>" or (len(ref_allele) > 1 and len(alt_allele) == 1):
        return "deletion"
    elif alt_allele == "<INS>" or (len(alt_allele) > 1 and len(ref_allele) == 1):
        return "insertion"
    elif alt_allele == "<INV>":
        return "inversion"
    elif alt_allele == "<DUP>":
        return "duplication"
    elif alt_allele == "<TDUP>":
        return "tandem duplication"
    elif alt_allele == "<CNV>":
        return "copy number variation"
    elif alt_allele.startswith("<"):
        # catch other special alleles
        return "complex"
    elif len(alt_allele) == 1 and len(ref_allele) == 1:
        return "snv"
    elif len(alt_allele) == len(ref_allele):
        return "mnv"
    elif len(alt_allele) > 1 and len(ref_allele) > 1:
        return "replacement"
    else:
        return "breakend"


def order_impact(df):
    order_impact = ["MODIFIER", "LOW", "MODERATE", "HIGH"]
    df["impact"] = pd.Categorical(df["impact"], order_impact)


def sort_calls(df):
    df.sort_values(snakemake.params.sorting, ascending=False, inplace=True)


def reorder_prob_cols(df):
    prob_df = df.filter(regex="^prob: ")
    if not prob_df.empty:
        prob_sums = prob_df.sum().sort_values(ascending=False)
        ordered_columns = prob_sums.keys()
        indexes = [df.columns.get_loc(col) for col in ordered_columns]
        start_index = min(indexes)
        updated_columns = df.columns.drop(labels=ordered_columns)
        for i, col in enumerate(ordered_columns):
            updated_columns = updated_columns.insert(start_index + i, col)
        df = df.reindex(columns=updated_columns)
    return df


def cleanup_dataframe(df):
    df = drop_low_prob_cols(df)
    df = reorder_prob_cols(df)
    df = format_floats(df)
    return df


calls = pd.read_csv(snakemake.input[0], sep="\t")
calls["clinical significance"] = (
    calls["clinical significance"]
    .apply(eval)
    .apply(sorted)
    .apply(", ".join)
    .replace("", np.nan)
)


if not calls.empty:
    # these below only work on non empty dataframes
    calls["vartype"] = calls.apply(get_vartype, axis="columns")
    order_impact(calls)
    sort_calls(calls)
else:
    calls["vartype"] = []


calls.set_index("gene", inplace=True, drop=False)
samples = get_samples(calls)


coding = ~pd.isna(calls["hgvsp"])
canonical = calls["canonical"]

noncoding_calls = calls[~coding & canonical]
noncoding_calls = cleanup_dataframe(noncoding_calls)
write(noncoding_calls, snakemake.output.noncoding)

coding_calls = calls[coding & canonical].drop(
    [
        "chromosome",
        "position",
        "reference allele",
        "alternative allele",
        "end position",
        "event",
        "id",
    ],
    axis="columns",
)
coding_calls = cleanup_dataframe(coding_calls)

# coding variants
# Here we have all variant info in hgvsp,
# hence we can drop the other variant specific columns.
write(
    coding_calls,
    snakemake.output.coding,
)

# TODO add possibility to also see non-canonical transcripts (low priority, once everything else works).
