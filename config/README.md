To configure this workflow, modify ``config/config.yaml`` according to your needs, following the explanations provided in the file.

Moreover, add your samples to the sample sheet ``config/samples.tsv`` and the corresponding sequencing units to ``config/units.tsv``.
All columns in the template sheets except ``purity`` are mandatory.

Finally, modify ``config/scenario.yaml`` according to your needs.
It allows to you to configure this pipeline for any arbitrary variant calling scenario (see the [Varlociraptor docs](https://varlociraptor.github.io/docs/calling/#generic-variant-calling)).
The file is templated with [Jinja2](https://jinja.palletsprojects.com) and will be rendered for each group of samples (see group column in the sample sheet).
For details, see the explanations in the file.