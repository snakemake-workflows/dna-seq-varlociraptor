import json
import re
import sys
import numpy as np
import pandas as pd

sys.stderr = open(snakemake.log[0], "w")


def write(df, path):
    df.drop(["canonical"], axis="columns", errors="ignore").dropna(how="all", axis="columns").to_csv(
        path, index=False, sep="\t"
    )


def format_floats(df):
    for col_name in df:
        if issubclass(df[col_name].dtype.type, np.floating):
            df[col_name] = [
                "{:.2e}".format(x) if x < 0.1 and x > 0 else round(x, 2)
                for x in df[col_name]
            ]
    return df


def trim_hgvsp_entries(df):
    if "hgvsp" in df.columns:
        df["hgvsp"] = df["hgvsp"].apply(
            lambda x: re.sub(r"^.*?:", "", x) if pd.notnull(x) else x
        )
    return df


def get_vaf_columns(df):
    return df.columns[df.columns.str.endswith(": allele frequency")]


def get_samples(df):
    return list(get_vaf_columns(df).str.replace(": allele frequency", ""))


def plot_data(df):
    samples = get_samples(df)
    vafcols = list(get_vaf_columns(df))
    varcol = "hgvsp"
    if "hgvsp" not in df.columns:
        varcol = "hgvsg"
    df = df[vafcols + [varcol]]
    df.columns = ["variant", *samples]
    return df


def plot_spec(samples):
    plots = [
        {
            "layer": [
                {"mark": "line"},
                {"mark": {"type": "point", "filled": True}},
                {
                    "mark": {"type": "text", "dx": 13, "align": "left"},
                    "encoding": {"text": {"field": "variant", "type": "nominal"}},
                    "transform": [{"filter": f"datum.sample == {samples[-1]}"}],
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
calls = format_floats(calls)
calls = trim_hgvsp_entries(calls)
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
