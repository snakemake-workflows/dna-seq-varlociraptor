import altair as alt
import pandas as pd

sys.stderr = open(snakemake.log[0], "w")

df = pd.read_csv(snakemake.input, sep="\t")
# Rename sample column into "frequencies"

chart = (
    alt.Chart(source)
    .mark_area()
    .encode(
        x="max_vaf:T",
        y=f"{snakemake.wildcards.group}:Q",
    )
)

chart.save(snakemake.output[0])
