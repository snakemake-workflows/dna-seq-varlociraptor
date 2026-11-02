import os
import sys

import polars as pl

sys.stderr = open(snakemake.log[0], "w")

species = snakemake.params.species
name = species.replace("_", " ").strip().lower()

names = pl.read_csv(
    snakemake.input.names,
    separator="|",
    has_header=False,
    new_columns=["taxid", "name", "unique_name", "name_class", "trailing"],
    schema_overrides={"taxid": pl.String, "name": pl.String, "name_class": pl.String},
    quote_char=None,
).with_columns(pl.col("taxid", "name", "name_class").str.strip_chars())

taxids = set(
    names.filter(
        (pl.col("name_class") == "scientific name")
        & (pl.col("name").str.to_lowercase() == name)
    )["taxid"].to_list()
)

if len(taxids) != 1:
    raise ValueError(
        f"Expected exactly one NCBI taxid for species '{species}', found: {sorted(taxids)}"
    )
taxid = taxids.pop()

# input file format of ganon build-custom: file, target, taxid
with open(snakemake.output[0], "w") as out:
    print(os.path.abspath(snakemake.input.genome), species, taxid, sep="\t", file=out)
