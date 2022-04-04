
# General settings
To configure this workflow, modify ``config/config.yaml`` according to your needs, following the explanations provided in the file.

# Sample sheet

Add samples to `config/samples.tsv`. For each sample, the columns `sample_name`, `alias`, `platform`, and `group` have to be defined. 
* Samples within the same `group` will be called jointly. 
* Aliases represent the name of the sample within its group (they can be the same as the sample name, or something simpler, e.g. tumor or normal).
* The `platform` column needs to contain the used sequencing plaform (one of 'CAPILLARY', 'LS454', 'ILLUMINA', 'SOLID', 'HELICOS', 'IONTORRENT', 'ONT', 'PACBIO’).

If mutational burdens shall be estimated for a sample, the to be used ``events`` from the calling scenario (see below) have to be specified in an additional column ``mutational_burden_events``. Multiple events have to be separated by commas within that column.

Missing values can be specified by empty columns or by writing `NA`. Lines can be commented out with `#`.

# Unit sheet

For each sample, add one or more sequencing units (runs, lanes or replicates) to the unit sheet `config/units.tsv`.
* Each unit has a `unit_name`, which can be e.g. a running number, or an actual run, lane or replicate id.
* Each unit has a `sample_name`, which associates it with the biological sample it comes from.
* For each unit, define either one (column `fq1`) or two (columns `fq1`, `fq2`) FASTQ files (these can point to anywhere in your system). 
* Alternatively, you can define an SRA (sequence read archive) accession (starting with e.g. ERR or SRR) by using a column `sra`. In the latter case, the pipeline will automatically download the corresponding paired end reads from SRA. If both local files and SRA accession are available, the local files will be preferred.
* Define adapters in the `adapters` column, by putting [cutadapt arguments](https://cutadapt.readthedocs.org) in quotation marks (e.g. `"-a ACGCGATCG -A GCTAGCGTACT"`).

Missing values can be specified by empty columns or by writing `NA`. Lines can be commented out with `#`.

# Calling scenario

Varlociraptor supports integrated uncertainty aware calling and filtering of variants for arbitrary scenarios. These are defined as so-called scenarios, via a [variant calling grammar](https://varlociraptor.github.io/docs/calling#generic-variant-calling).
* For each group, a scenario is rendered via [Jinja](https://jinja.palletsprojects.com).
* Therefore, edit the template scenario (`scenario.yaml`) according to your needs. The sample sheet is available for jinja rendering as a pandas data frame in the variable `samples`. This allows to customize the scenario according to the contents of the sample sheet. You can therefore add additional columns to the sample sheet (e.g. purity) and access them in the scenario template, in order to pass the information to Varlociraptor.

# Primer trimming

For panel data the pipeline allows trimming of amplicon primers on both ends of a fragment but also on a single end only. 
In case of single end primers these are supposed to be located at the left end of a read.
When primer trimming is enabled, primers have to be defined either directly in the `config.yaml` or in a seperate tsv-file.
Defining primers directly in the config file is prefered when all samples come from the same primer set.
In case of different panels, primers have to be set panel-wise in a seperate tsv-file.
For each panel the following columns need to be set: `panel`, `fa1` and `fa2` (optional).
Additionally, for each sample the corresponding panel must be defined in `samples.tsv` (column `panel`).
For single primer trimming only, the first entry in the config (respective in the tsv file) needs to be defined.
