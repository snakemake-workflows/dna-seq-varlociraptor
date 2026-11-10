from pathlib import Path

import polars as pl
import pysam


def get_colname(col):
    if col:
        return str(Path(col).stem)
    else:
        return "species"


def get_num_reads(stats):
    # extract total number of reads from bam stats file
    with open(stats) as f:
        for line in f:
            if line.startswith("SN\tsequences:"):
                return int(line.strip().split()[-1])
    raise ValueError(f"Could not find number of sequences in {stats}")


num_reads = pl.DataFrame({"sample": snakemake.params.samples, "num_reads": [get_num_reads(stats) for stats in snakemake.input.stats]})

data = pl.read_csv(snakemake.input[0], separator="\t")
data.rename({col: get_colname(col) for col in data.columns}).transpose(include_header=True).join(num_reads, on="sample").select(
    pl.col("species"),
    pl.exclude(["species", "num_reads"]) / pl.col("num_reads"),
).write_parquet(snakemake.output[0])