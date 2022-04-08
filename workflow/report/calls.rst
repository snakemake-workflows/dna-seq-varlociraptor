{% set event = snakemake.config['calling']['fdr-control']['events'][snakemake.wildcards.event] %}
{% if "desc" in event %}
{{ event["desc"] }}
{% else %}
Variants for the **event {{snakemake.wildcards.event}}**.
{% endif %}

{% if snakemake.config["calling"]["fdr-control"].get("local") %}
{% set fdr = "`**local** false discovery rate <https://en.wikipedia.org/wiki/False_discovery_rate>`_, i.e. the probability for each variant to be a false positive under consideration of above event definition is at most" %}
{% else %}
{% set fdr = "`false discovery rate <https://en.wikipedia.org/wiki/False_discovery_rate>`_, i.e. the expected fraction of false positives according to the given event definition among all variants in this callset is at most" %}
{% endif %}

The {{ fdr }} {{snakemake.config["calling"]["fdr-control"]["threshold"]}}, based on the posterior probabilities delivered by `Varlociraptor's <https://varlociraptor.github.io>`_ Bayesian latent variable model.

The calculated probabilties entail various sources of uncertainty like typing and mapping uncertainty, as well as typical biases like strand bias, read orientation bias, read position bias, homopolymer errors, and more.
They are further capable of properly capturing the uncertainty increase occuring with lower read depths or allele frequencies, so that no read depth or minimum allele frequency filtering beyond the Varlociraptor scenario configuration is necessary.