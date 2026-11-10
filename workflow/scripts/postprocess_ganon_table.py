from pathlib import Path

import polars as pl


def get_colname(col):
    if col:
        return str(Path(col).stem)
    else:
        return "species"


data = pl.read_csv(snakemake.input[0], separator="\t")
data.rename({col: get_colname(col) for col in data.columns}).write_parquet(snakemake.output[0])