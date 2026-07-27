import polars as pl
import networkx as nx
import altair as alt

samples = pl.from_pandas(snakemake.params.samples).select("sample_name", "group")

# read pair data and restrict to pairs that have a high concordance
pairs = (
    pl.read_csv(snakemake.input[0], separator="\t")
    .rename({"#sample_a": "sample_a"})
    .join(samples, left_on="sample_a", right_on="sample_name", how="left", suffix="_a")
    .join(samples, left_on="sample_b", right_on="sample_name", how="left", suffix="_b")
    .filter(
        (pl.col("relatedness") >= 0.9)
        | (pl.col("concordance") >= 0.9)
        | (pl.col("group") == pl.col("group_b"))
    )
)

# obtain graph
graph = nx.from_pandas_edgelist(pairs, source="sample_a", target="sample_b")

# reduce to connected components that are impure
impure_subset = set()
for component in nx.connected_components(graph):
    if (
        samples.filter(pl.col("sample_name").is_in(component))
        .get_column("group")
        .unique()
        .len()
        > 1
    ):
        impure_subset.update(component)

# generate layout
layout = nx.spring_layout(graph.subgraph(impure_subset), seed=2798791)
coords = pl.DataFrame(
    {
        "sample_name": layout.keys(),
        "x": [pos[0] for pos in layout.values()],
        "y": [pos[1] for pos in layout.values()],
    }
)

# add coords to sample sheet and only keep samples that are of interest
data = samples.join(coords, how="inner", on="sample_name")

# add edges for visualization and encode strength
edges = (
    pairs.filter(
        pl.col("sample_a").is_in(impure_subset) & pl.col("sample_b").is_in(impure_subset)
    )
    .select(
        "sample_a",
        "sample_b",
        "relatedness",
        "concordance",
        pl.max_horizontal("relatedness", "concordance").alias("edge_strength"),
    )
    .join(
        coords.rename({"sample_name": "sample_a", "x": "x_a", "y": "y_a"}),
        on="sample_a",
        how="inner",
    )
    .join(
        coords.rename({"sample_name": "sample_b", "x": "x_b", "y": "y_b"}),
        on="sample_b",
        how="inner",
    )
)

# plot
base = alt.Chart(data).encode(alt.X("x").axis(None), alt.Y("y").axis(None))
edge_base = alt.Chart(edges).encode(
    alt.X("x_a").axis(None),
    alt.Y("y_a").axis(None),
    alt.X2("x_b"),
    alt.Y2("y_b"),
)
(
    edge_base.mark_line(clip=False, color="black")
    .encode(
        alt.Opacity("edge_strength").scale(domain=[0, 1]),
        tooltip=[
            alt.Tooltip("sample_a"),
            alt.Tooltip("sample_b"),
            alt.Tooltip("relatedness"),
            alt.Tooltip("concordance"),
            alt.Tooltip("edge_strength"),
        ],
    )
    + base.mark_circle(tooltip=True, clip=False, size=30).encode(
        alt.Color("group").scale(scheme="category20")
    )
    + base.mark_text(dx=5, dy=5, clip=False, align="left").encode(
        alt.Text("sample_name")
    )
).interactive().properties(
    width="container",
    height=700,
).configure_view(
    stroke=None
).save(
    snakemake.output[0]
)
