import polars as pl
import altair as alt
alt.data_transformers.enable("vegafusion")


data = pl.read_parquet(snakemake.input[0]).unpivot(index="sample", variable_name="species", value_name="fraction").filter(pl.col("fraction") >= 0.001)

alt.Chart(data).mark_circle(tooltip=True).encode(
    alt.X("sample:N"),
    alt.Y("fraction:Q"),
    alt.Color("species:N")
).interactive().save(snakemake.output[0])