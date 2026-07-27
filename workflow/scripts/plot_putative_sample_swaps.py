import polars as pl
import networkx as nx
import altair as alt

samples = pl.from_pandas(snakemake.params.samples).select("sample_name", "group")

pairs = (
    pl.read_csv(snakemake.input[0], separator="\t")
    .rename({"#sample_a": "sample_a"})
    .join(samples, left_on="sample_a", right_on="sample_name", how="left", suffix="_a")
    .join(samples, left_on="sample_b", right_on="sample_name", how="left", suffix="_b")
    .filter(
        (pl.col("relatedness") >= 0.9)
        & (pl.col("concordance") >= 0.9)
        & (pl.col("group") == pl.col("group_b"))
    )
)

graph = nx.from_pandas_edgelist(pairs, source="sample_a", target="sample_b")

layout = nx.spring_layout(graph, seed=2798791)

coords = pl.DataFrame(
    {
        "sample_name": layout.keys(),
        "x": [pos[0] for pos in layout.values()],
        "y": [pos[1] for pos in layout.values()],
    }
).cast({"sample_name": str})
breakpoint()

base = alt.Chart(
    samples.join(coords, on="sample_name")
).encode(alt.X("x"), alt.Y("y"))

(
    base.mark_circle(tooltip=True).encode(
        alt.Color("group")
    ) + base.mark_text(dx=5, dy=5).encode(alt.Text("sample"))
).save(snakemake.output[0])
