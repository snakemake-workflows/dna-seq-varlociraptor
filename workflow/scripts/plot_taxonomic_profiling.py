import polars as pl
import altair as alt

alt.data_transformers.enable("vegafusion")


def top_n_species(data, n=19):
    top_species = (
        data.group_by("species")
        .agg(pl.max("fraction").alias("max_fraction"))
        .sort("max_fraction", descending=True)
        .head(n)
        .get_column("species")
    )
    return top_species


def plot_organism_group(data, organism_group):
    data = data.filter(pl.col("organism_group") == organism_group)
    data = data.filter(pl.col("species").is_in(top_n_species(data, 100))).with_columns(
        pl.when(pl.col("species").is_in(top_n_species(data, 19)))
        .then(pl.col("species"))
        .otherwise(pl.lit("Other"))
        .alias("top_species")
    )

    species_selection = alt.selection_point(fields=["top_species"], bind="legend")

    return (
        alt.Chart(data)
        .mark_circle(tooltip=True, size=50)
        .encode(
            alt.X("sample:N"),
            alt.Y("fraction:Q"),
            alt.Color("top_species:N").title(organism_group).scale(scheme="category20"),
            # in the tooltip, show the actual species, not just the top_species category
            alt.Tooltip(["sample:N", "species:N", "fraction:Q", "count:Q"]),
            opacity=alt.when(species_selection)
            .then(alt.value(1))
            .otherwise(alt.value(0.1)),
        )
        .add_params(species_selection)
        .interactive()
    )


data = pl.read_parquet(snakemake.input[0]).filter(pl.col("species") != "unclassified")

alt.vconcat(
    *(
        plot_organism_group(data, organism_group)
        for organism_group in data.get_column("organism_group").unique()
    )
).resolve_scale(color="independent", x="shared").save(snakemake.output[0])
