import sys
from jinja2 import Template
import pandas as pd
import fastparquet

sys.stderr = open(snakemake.log[0], "w")


with open(snakemake.input[0]) as template, open(snakemake.output[0], "w") as out:
    samples = pd.read_parquet(snakemake.params.samples)
    out.write(
        Template(template.read()).render(
            samples=samples[samples["group"] == snakemake.wildcards.group]
        )
    )
