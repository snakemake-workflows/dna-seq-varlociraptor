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


# Drop columns, where both breakpoints will always have identical annotation,
# because varlociraptor evaluates them together as one event.
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
    snakemake.input["varlociraptor"],
    sep="\t",
    dtype={
        "exon": "Int64",
        "chromosome": str,
        "position": str,
    },
)
calls[["feature_name", "feature_id"]] = calls[["feature_id", "feature_name"]].fillna(
    "('not in exon',)"
)
calls = calls.drop(["reference allele", "alternative allele"], axis="columns")
# Obtain first entry of columns annotated as tuple by arriba
for col in ["feature_id", "feature_name"]:
    calls[col] = calls[col].apply(eval).apply(lambda x: x[0])

first_partners = calls["id"].str.contains("a$")
first_calls = calls[first_partners]
# The dropped columns are identical between breakpoints of the same event.
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

# From this file, we bring back useful arriba annotations for the final
# datavzrd table for the report.
annotated_candidates = pd.read_csv(
    snakemake.input["arriba"][0],
    sep="\t",
    usecols=[
        "gene_id1",
        "breakpoint1",
        "transcript_id1",
        "site1",
        "gene_id2",
        "breakpoint2",
        "transcript_id2",
        "site2",
        "type",
        "reading_frame",
    ],
    dtype=str
)
# For matching entries on chromosome and position, we need to split this column.
annotated_candidates[["chrom1", "pos1"]] = annotated_candidates["breakpoint1"].str.split(":", expand=True)
annotated_candidates[["chrom2", "pos2"]] = annotated_candidates["breakpoint2"].str.split(":", expand=True)
annotated_candidates = annotated_candidates.drop(columns=["breakpoint1", "breakpoint2"])

# We merge only on chromosome and position, as some breakpoints don't have gene
# name annotations and will receive a 'not in exon' entry (see above). Thus,
# their gene designation will differ from those provided by the arriba
# annotation.
paired_fusions_with_arriba_annotations = paired_fusions.merge(
    right=annotated_candidates,
    how="left",
    left_on=[
        "chromosome_partner1",
        "position_partner1",
        "chromosome_partner2",
        "position_partner2",
    ],
    right_on=[
        "chrom1",
        "pos1",
        "chrom2",
        "pos2",
    ],
)
# Some breakpoints have multiple gene annotations, and in merging only by
# chromosome and position this will lead to entries that are the cartesian
# product of the entries in the two sets. While keeping entries with a
# 'not in exon` feature_name, we filter those duplicates down to entries
# where the feature_name and gene_id agree.
paired_fusions_with_arriba_annotations = paired_fusions_with_arriba_annotations.loc[
    (paired_fusions_with_arriba_annotations['feature_name_partner1'] == 'not in exon') |
    (paired_fusions_with_arriba_annotations['feature_name_partner2'] == 'not in exon') |
    (
        ( paired_fusions_with_arriba_annotations['feature_name_partner1'] == paired_fusions_with_arriba_annotations['gene_id1']) &
        ( paired_fusions_with_arriba_annotations['feature_name_partner2'] == paired_fusions_with_arriba_annotations['gene_id2'])
    )
# Remove duplicate columns only needed for record matching.
].drop(
    columns=[
        "gene_id1",
        "chrom1",
        "pos1",
        "gene_id2",
        "chrom2",
        "pos2",
    ]
)

write(paired_fusions_with_arriba_annotations, snakemake.output.fusions)
