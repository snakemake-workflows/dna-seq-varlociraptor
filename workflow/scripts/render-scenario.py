from jinja2 import Template
import pandas as pd

with open(snakemake.input[0]) as template, open(snakemake.output[0], "w") as out:
    samples = snakemake.params.samples
    out.write(Template(template.read()).render(
        samples=samples[samples["group"] == snakemake.wildcards.group]
    ))
