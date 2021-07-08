Estimated mutational burden for sample {{ snakemake.wildcards.sample }}, 
using events {{ snakemake.params.events.split(" ")|join(", ") }}.

The calculated plot shows the expected number of mutations (derived 
from the posterior event probabilities calculated by Varlociraptor) 
versus minimum variant allele frequencies.