import altair as alt
import pandas as pd

sys.stderr = open(snakemake.log[0], "w")

df = pd.read_csv(snakemake.input[0], sep="\t")

chart = (
    alt.Chart(df)
    .mark_area(interpolate="basis")
    .encode(
        x=alt.X("min_vaf:Q", scale=alt.Scale(reverse=True)),
        y=f"Frequency:Q",
        color="Signature:N"
    )
)

chart.save(snakemake.output[0])
