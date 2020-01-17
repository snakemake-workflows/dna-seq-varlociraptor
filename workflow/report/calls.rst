Variants for the **event {{snakemake.wildcards.event}}**, filtered by controlling false discovery rate at {{snakemake.config["calling"]["fdr-control"]["threshold"]}}.

The top left panel shows primary information about each variant.
Upon click on the variant, the right panel shows details about the impact, as well as all levels of evidence from event probabilities, to posterior allele frequency distributions down to likelihoods at read/fragment level. The latter is shown as `Kass/Raftery-scored Bayes Factors <https://en.wikipedia.org/wiki/Bayes_factor>`_, denoting the strength of evidence for the variant allele.
