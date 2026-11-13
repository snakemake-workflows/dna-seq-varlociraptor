import polars as pl
import altair as alt

alt.data_transformers.enable("vegafusion")


data = (
    pl.read_parquet(snakemake.input[0])
    .unpivot(index="sample", variable_name="species", value_name="fraction")
    .filter((pl.col("fraction") >= 0.0001) & (pl.col("species") != "unclassified"))
)

# restrict to top 20 and color all the rest as other
top_species = (
    data.group_by("species")
    .agg(pl.max("fraction").alias("max_fraction"))
    .sort("max_fraction", descending=True)
    .head(19)
    .get_column("species")
)
print(data.filter(pl.col("species").is_in(top_species)))
data = data.with_columns(
    pl.when(pl.col("species").is_in(top_species))
    .then(pl.col("species"))
    .otherwise(pl.lit("Other"))
    .alias("top_species")
)

alt.Chart(data).mark_circle(tooltip=True).encode(
    alt.X("sample:N"),
    alt.Y("fraction:Q"),
    alt.Color("top_species:N").title("species"),
    # in the tooltip, show the actual species, not just the top_species category
    alt.Tooltip(["sample:N", "species:N", "fraction:Q"]),
).interactive().save(snakemake.output[0])
