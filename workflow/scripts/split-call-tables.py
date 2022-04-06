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
        df = df.dropna(how="all", axis="columns")
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


def drop_zero_vaf_cols(df):
    return drop_cols_by_predicate(
        df,
        df.columns[df.columns.str.endswith(": allele frequency")],
        lambda vafs: vafs == 0.0,
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
    else:
        return "breakend"


def plot_data(df):
    samples = get_samples(df)
    vafcols = list(get_vaf_columns(df))
    varcol = "hgvsp"
    if pd.isna(df["hgvsp"]).all():
        varcol = "hgvsg"
    df = df[vafcols + [varcol]]
    df.columns = [*samples, "variant"]
    return df


def sort_calls(df):
    order_impact = ["MODIFIER", "LOW", "MODERATE", "HIGH"]
    df["impact"] = pd.Categorical(df["impact"], order_impact)
    df.sort_values(snakemake.params.sorting, ascending=False, inplace=True)
    return df


def plot_spec(samples):
    plot_spec = {
        "$schema": "https://vega.github.io/schema/vega-lite/v5.json",
        "description": "Drag the sliders to highlight points.",
        "transform": [{"fold": samples, "as": ["sample", "vaf"]}],
        "mark": "rect",
        "encoding": {
            "x": {"field": "sample", "type": "ordinal"},
            "color": {
                "field": "vaf",
                "type": "quantitative",
                "scale": {"domain": [0, 1]},
            },
            "y": {"field": "variant"},
            "href": {"field": "variant-link", "type": "nominal"},
        },
    }
    return plot_spec


def plot_spec_old(samples):
    plots = [
        {
            "layer": [
                {"mark": "line"},
                {"mark": {"type": "point", "filled": True}},
                {
                    "mark": {"type": "text", "dx": 13, "align": "left"},
                    "encoding": {"text": {"field": "variant", "type": "nominal"}},
                    "transform": [{"filter": f"datum.sample == '{samples[-1]}'"}],
                },
            ],
            "encoding": {
                "x": {"field": "sample", "type": "ordinal", "sort": samples},
                "y": {
                    "field": "vaf",
                    "type": "quantitative",
                    "scale": {"domain": [0, 1.0]},
                },
                "color": {"field": "variant", "legend": None},
                "href": {"field": "variant-link"},
                "opacity": {
                    "condition": {
                        "test": {
                            "and": [{"param": f"{sample}brush"} for sample in samples]
                        },
                        "value": 1.0,
                    },
                    "value": 0.2,
                },
            },
        }
    ]

    for sample in samples:
        plots.append(
            {
                "mark": "tick",
                "encoding": {
                    "x": {
                        "field": sample,
                        "type": "quantitative",
                        "scale": {"domain": [0, 1.0]},
                    }
                },
                "params": [
                    {
                        "name": f"{sample}brush",
                        "select": {"type": "interval", "encodings": ["x"]},
                    }
                ],
            }
        )

    plot_spec = {
        "$schema": "https://vega.github.io/schema/vega-lite/v5.json",
        "description": "Drag the sliders to highlight points.",
        "transform": [{"fold": samples, "as": ["sample", "vaf"]}],
        "vconcat": plots,
    }
    return plot_spec


calls = pd.read_csv(snakemake.input[0], sep="\t")
calls["clinical significance"] = (
    calls["clinical significance"]
    .apply(eval)
    .apply(sorted)
    .apply(", ".join)
    .replace("", np.nan)
)
calls["vartype"] = calls.apply(get_vartype, axis="columns")
calls = sort_calls(calls)
calls = format_floats(calls)
calls = drop_low_prob_cols(calls)
calls = drop_zero_vaf_cols(calls)
calls.set_index("gene", inplace=True, drop=False)
samples = get_samples(calls)


coding = ~pd.isna(calls["hgvsp"])
canonical = calls["canonical"]

noncoding_calls = calls[~coding & canonical]

coding_calls = calls[coding & canonical].drop(
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
)


write(noncoding_calls, snakemake.output.noncoding)

# coding variants
# Here we have all variant info in hgvsp,
# hence we can drop the other variant specific columns.
write(
    coding_calls,
    snakemake.output.coding,
)


write(plot_data(coding_calls), snakemake.output.coding_plot_data)
write(plot_data(noncoding_calls), snakemake.output.noncoding_plot_data)

with open(snakemake.output.plot_spec, "w") as out:
    json.dump(plot_spec(samples), out)

# TODO add possibility to also see non-canonical transcripts (low priority, once everything else works).
