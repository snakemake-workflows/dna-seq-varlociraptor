import sys

sys.stderr = open(snakemake.log[0], "w")

import pandas as pd


def write(df, path):
    if not df.empty:
        remaining_columns = df.dropna(how="all", axis="columns").columns.tolist()
        df = df[remaining_columns]
    df.to_csv(path, index=False, sep="\t")


def get_prefix_columns(df, prefix):
    return df.columns[df.columns.str.startswith(prefix)]


def get_samples(df):
    return list(
        df.columns[df.columns.str.endswith(": allele frequency")].str.replace(
            ": allele frequency", ""
        )
    )


def drop_shared_columns(df):
    prob_columns = get_prefix_columns(df, "prob: ")
    df = df.drop(prob_columns, axis="columns")
    samples = get_samples(df)
    sample_columns = [
        column for sample in samples for column in get_prefix_columns(df, sample)
    ]
    return df.drop(sample_columns, axis="columns")


def join_short_obs(df, samples):
    for sample in samples:
        sobs_loc = df.columns.get_loc(f"{sample}: observations")
        df.insert(
            sobs_loc,
            f"{sample}: short observations",
            df[f"{sample}: short ref observations"]
            + ","
            + df[f"{sample}: short alt observations"],
        )
    df = df.drop(
        df.columns[
            df.columns.str.endswith(": short ref observations")
            | df.columns.str.endswith(": short alt observations")
        ],
        axis=1,
    )
    return df


calls = pd.read_csv(
    snakemake.input[0],
    sep="\t",
    dtype={"exon": "Int64"},
)
calls[["feature_name", "feature_id"]] = calls[["feature_id", "feature_name"]].fillna(
    "('intronic',)"
)
calls = calls.drop(["reference allele", "alternative allele"], axis="columns")
# Obtain first entry of columns annotated as tuple by arriba
for col in ["feature_id", "feature_name"]:
    calls[col] = calls[col].apply(eval).apply(lambda x: x[0])

first_partners = calls["id"].str.contains("a$")
first_calls = calls[first_partners]
first_calls = drop_shared_columns(first_calls)
second_calls = calls[~first_partners]
samples = get_samples(second_calls)
second_calls = join_short_obs(second_calls, samples)
paired_fusions = first_calls.merge(
    second_calls,
    left_on="mateid",
    right_on="id",
    suffixes=("_partner1", "_partner2"),
)
paired_fusions = paired_fusions.filter(regex="^(?!mateid|id)")

write(paired_fusions, snakemake.output.fusions)
