import altair as alt
import pandas as pd

sys.stderr = open(snakemake.log[0], "w")

signatures_df = pd.read_csv(snakemake.input.signatures, sep="\t")

signatures = (
    alt.Chart(signatures_df)
    .mark_area(interpolate="monotone")
    .encode(
        x=alt.X("Minimum VAF:Q", scale=alt.Scale(reverse=True)),
        y=f"Frequency:Q",
        color="Signature:N",
        tooltip="Description:N",
    )
)

mut_counts_df = pd.read_csv(snakemake.input.counts, sep="\t")

counts = (
    alt.Chart(mut_counts_df)
    .mark_line(interpolate="basis", color="black")
    .encode(
        x=alt.X("Minimum VAF:Q", scale=alt.Scale(reverse=True)),
        y="Mutation Count:Q",
    )
)

final_chart = alt.layer(signatures, counts).resolve_scale(y="independent")

final_chart.save(snakemake.output[0])
