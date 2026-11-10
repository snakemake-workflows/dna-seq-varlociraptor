import polars as pl
import altair as alt
alt.data_transformers.enable("vegafusion")


data = pl.read_parquet(snakemake.input[0]).unpivot(index="sample", variable_name="species", value_name="fraction").filter((pl.col("fraction") >= 0.001) & (pl.col("species") != "Homo sapiens"))

# restrict to top 20 and color all the rest as other
top_species = data.groupby("species").agg(pl.sum("fraction").alias("total_fraction")).sort("total_fraction", reverse=True).head(20)["species"]
data = data.with_column(
    pl.when(pl.col("species").is_in(top_species))
    .then(pl.col("species"))
    .otherwise("Other")
    .alias("top_species")
)

alt.Chart(data).mark_circle(tooltip=True).encode(
    alt.X("sample:N"),
    alt.Y("fraction:Q"),
    alt.Color("top_species:N").title("species"),
    # in the tooltip, show the actual species, not just the top_species category
    alt.Tooltip(["sample:N", "species:N", "fraction:Q"]),
).interactive().save(snakemake.output[0])