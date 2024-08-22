import sys

sys.stderr = open(snakemake.log[0], "w")

import numpy as np
import pandas as pd
import pysam

from typing import Generator

PROB_EPSILON = 0.01  # columns with all probabilities below will be dropped


def write(df, path):
    df["mane_plus_clinical"][df["mane_plus_clinical"].notna()] = True
    if not df.empty:
        remaining_columns = df.dropna(how="all", axis="columns").columns.tolist()
        if path == snakemake.output.coding:
            # ensure that these columns are kept, even if they contain only NAs in a coding setting
            remaining_columns.extend(["revel", "hgvsp", "symbol"])
            remaining_columns = [col for col in df.columns if col in remaining_columns]
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


def reorder_vaf_cols(df):
    vaf_columns = df.filter(regex=": allele frequency$").columns
    other_columns = df.drop(columns=vaf_columns).columns
    split_index = other_columns.get_loc("consequence")
    left_columns = other_columns[0:split_index]
    right_columns = other_columns[split_index:]
    reordered_columns = left_columns.append([vaf_columns, right_columns])
    return df[reordered_columns]


def cleanup_dataframe(df):
    df = drop_low_prob_cols(df)
    df = reorder_prob_cols(df)
    df = reorder_vaf_cols(df)
    df = format_floats(df)
    return df


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


def bin_max_vaf(df, samples):
    af_columns = [f"{sample}: allele frequency" for sample in samples]
    max_vaf = df[af_columns].apply("max", axis=1)
    df["binned max vaf"] = pd.cut(
        max_vaf, [0, 0.33, 0.66, 1.0], labels=["low", "medium", "high"]
    )
    df["binned max vaf"] = pd.Categorical(
        df["binned max vaf"], ["low", "medium", "high"]
    )
    return df


class PopulationDb:
    def __init__(self, path):
        self.contig = None
        self.pos = None
        self._variants = None
        self.bcf = pysam.VariantFile(path)

    def variants(
        self, contig: str, pos: int, alt: str
    ) -> Generator[pysam.VariantRecord, None, None]:
        """Return variants at given position"""
        if not self._is_in_interval(contig, pos):
            self.contig = contig
            self.pos = pos
            self._variants = self._load_variants()
        for variant in self._variants:
            if variant.pos == pos and variant.alts[0] == alt:
                yield variant
            if variant.pos > pos:
                break

    def annotate_row(self, row):
        # TODO: deal with SVs
        db_vars = self.variants(
            row["chromosome"], row["position"], row["alternative allele"]
        )
        return ",".join(
            [
                f"{name}:{sample['AF'][0]:0.2f}"
                for variant in db_vars
                for name, sample in zip(
                    self.bcf.header.samples, variant.samples.values()
                )
                if sample["AF"][0] and sample["AF"][0] > 0.0
            ]
        )

    def _load_variants(self):
        return self.bcf.fetch(str(self.contig), self.pos, self.end)

    @property
    def end(self):
        return self.pos + 1000

    def _is_in_interval(self, contig: str, pos: int):
        return (
            self.pos is not None
            and self.contig == contig
            and self.pos <= pos <= self.end
        )


def select_spliceai_effect(calls):
    spliceai_columns = calls.filter(like="spliceai", axis=1)
    max_spliceai_effects = (
        spliceai_columns.idxmax(axis=1).astype(str).str.removeprefix("spliceai ") + ": "
    )
    max_score = spliceai_columns.max(axis=1)
    # Do not annotate effect for variants with high uncertainty (prob below 0.2)
    max_spliceai_effects[max_score < 0.2] = np.nan
    col_index = calls.columns.get_loc("spliceai acceptor gain")
    calls = calls.drop(calls.filter(like="spliceai", axis=1).columns, axis=1)
    calls.insert(col_index, "spliceai_effect", max_spliceai_effects)
    calls.insert(col_index, "spliceai", max_score)
    return calls


calls = pd.read_csv(snakemake.input[0], sep="\t")
calls["clinical significance"] = (
    calls["clinical significance"]
    .apply(eval)
    .apply(sorted)
    .apply(", ".join)
    .replace("", np.nan)
)

calls["protein alteration (short)"] = (
    calls["protein alteration (short)"].apply(eval).apply("/".join).replace("", np.nan)
)

samples = get_samples(calls)

if calls.columns.str.endswith(": allele frequency").any():
    calls = bin_max_vaf(calls, samples)

if snakemake.input.population_db and not calls.empty:
    population_db = PopulationDb(snakemake.input.population_db)
    calls["population"] = calls.apply(
        population_db.annotate_row, axis="columns"
    ).replace("", np.nan)

if not calls.empty:
    # these below only work on non empty dataframes
    calls["vartype"] = calls.apply(get_vartype, axis="columns")
    order_impact(calls)
    sort_calls(calls)
else:
    calls["vartype"] = []


calls.set_index("gene", inplace=True, drop=False)

if calls.columns.str.endswith(": short ref observations").any():
    calls = join_short_obs(calls, samples)

if calls.columns.str.startswith("spliceai").any():
    calls = select_spliceai_effect(calls)

coding = ~pd.isna(calls["hgvsp"])
canonical = calls["canonical"].notnull()
mane_plus_clinical = calls["mane_plus_clinical"].notnull()
canonical_mane = canonical | mane_plus_clinical

noncoding_calls = calls[~coding & canonical_mane]
noncoding_calls = cleanup_dataframe(noncoding_calls)
write(noncoding_calls, snakemake.output.noncoding)

coding_calls = calls[coding & canonical_mane].drop(
    [
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
# TODO add possibility to also see non-canoncical or non-mane+clinical transcripts (low priority, once everything else works).
