
# General settings
To configure this workflow, modify ``config/config.yaml`` according to your needs, following the explanations provided in the file.

# Sample and unit sheet

* Add samples to `config/samples.tsv`. For each sample, the columns `sample_name`, `alias`, `platform`, and `group` have to be defined. Samples within the same `group` will be called jointly. Aliases represent the name of the sample within its group (they can be the same as the sample name, or something simpler, e.g. tumor or normal).
* For each sample, add one or more sequencing units (runs, lanes or replicates) to the unit sheet `config/units.tsv`. For each unit, define adapters, and either one (column `fq1`) or two (columns `fq1`, `fq2`) FASTQ files (these can point to anywhere in your system). Alternatively, you can define an SRA (sequence read archive) accession (starting with e.g. ERR or SRR) by using a column `sra`. In the latter case, the pipeline will automatically download the corresponding paired end reads from SRA. If both local files and SRA accession are available, the local files will be preferred.

Missing values can be specified by empty columns or by writing `NA`.

# Calling scenario

Varlociraptor supports integrated uncertainty aware calling and filtering of variants for arbitrary scenarios. These are defined as so-called scenarios, via a [variant calling grammar](https://varlociraptor.github.io/docs/calling#generic-variant-calling).
* For each group, a scenario is rendered via [Jinja](https://jinja.palletsprojects.com).
* Therefore, edit the template scenario (`scenario.yaml`) according to your needs. The sample sheet is available for jinja rendering as a pandas data frame in the variable `samples`. This allows to customize the scenario according to the contents of the sample sheet. You can therefore add additional columns to the sample sheet (e.g. purity) and access them in the scenario template, in order to pass the information to Varlociraptor.
