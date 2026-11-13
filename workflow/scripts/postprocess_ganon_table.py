from pathlib import Path

import polars as pl
import pysam


def get_colname(col):
    if col:
        return Path(col).stem
    else:
        return "species"


def get_num_reads(stats):
    # extract total number of reads from bam stats file
    with open(stats) as f:
        for line in f:
            if line.startswith("SN\tsequences:"):
                return int(line.strip().split()[-1])
    raise ValueError(f"Could not find number of sequences in {stats}")


def get_organism_groups(taxpath, db_spec):
    if db_spec["category"] == "organism_groups":
        group = db_spec["spec"].split(":")[0]
    else:
        group = "custom"

    return pl.read_csv(
        taxpath,
        separator="\t",
        columns=[3],
        has_header=False,
        new_columns=["species"],
    ).with_columns(pl.lit(group).alias("organism_group"))


organism_groups = pl.concat(
    [
        get_organism_groups(taxpath, db_spec)
        for taxpath, db_spec in zip(snakemake.input.tax, snakemake.params.db_specs)
    ]
)
num_reads = pl.DataFrame(
    {
        "sample": snakemake.params.samples,
        "num_reads": [get_num_reads(stats) for stats in snakemake.input.stats],
    }
)

data = pl.read_csv(snakemake.input.profiling, separator="\t")
data.rename({col: get_colname(col) for col in data.columns}).transpose(
    column_names="species", include_header=True, header_name="sample"
).join(num_reads, on="sample").select(
    pl.col("sample"),
    pl.exclude(["sample", "num_reads"]) / pl.col("num_reads"),
).unpivot(index="sample", variable_name="species", value_name="fraction").join(
    organism_groups, on="species", how="left"
).write_parquet(snakemake.output[0])
