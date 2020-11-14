# Snakemake workflow: DNA-seq variant calling with Varlociraptor

[![Snakemake](https://img.shields.io/badge/snakemake-≥5.7.0-brightgreen.svg)](https://snakemake.bitbucket.io)
[![GitHub actions status](https://github.com/snakemake-workflows/dna-seq-varlociraptor/workflows/Tests/badge.svg?branch=master)](https://github.com/snakemake-workflows/dna-seq-varlociraptor/actions?query=branch%3Amaster+workflow%3ATests)

This workflow detects genomic variants with [Delly](https://github.com/dellytools/delly) and [Freebayes](https://github.com/ekg/freebayes), followed by statistical assessment with [Varlociraptor](https://varlociraptor.github.io). It is designed to flexibly define calling groups, and directly integrates the fetching of SRA samples (if required) and reference data (the latter making use of [between workflow caching](https://snakemake.readthedocs.io/en/stable/executing/caching.html)).

**Note:** at the moment, [Varlociraptor](https://varlociraptor.github.io) is limited to SNVs, MNVs, small and large (structural) indels and hence also this workflow. This will change with future releases of [Varlociraptor](https://varlociraptor.github.io).

## Authors

* Felix Mölder (@FelixMoelder)
* Johannes Köster (@johanneskoester)

## Usage

In any case, if you use this workflow in a paper, don't forget to give credits to the authors by citing the URL of this (original) repository.

#### Step 1: Obtain a copy of this workflow

1. Create a new github repository using this workflow [as a template](https://help.github.com/en/articles/creating-a-repository-from-a-template).
2. [Clone](https://help.github.com/en/articles/cloning-a-repository) the newly created repository to your local system, into the place where you want to perform the data analysis.

#### Step 2: Configure workflow

##### General settings
Configure the workflow according to your needs via editing the file `config.yaml`.

##### Sample and unit sheet

* Add samples to `config/samples.tsv`. For each sample, the columns `sample_name`, `alias`, `platform`, and `group` have to be defined. Samples within the same `group` will be called jointly. Aliases represent the name of the sample within its group (they can be the same as the sample name, or something simpler, e.g. tumor or normal).
* For each sample, add one or more sequencing units (runs, lanes or replicates) to the unit sheet `config/units.tsv`. For each unit, define adapters, and either one (column `fq1`) or two (columns `fq1`, `fq2`) FASTQ files (these can point to anywhere in your system). Alternatively, you can define an SRA (sequence read archive) accession (starting with e.g. ERR or SRR) by using a column `sra`. In the latter case, the pipeline will automatically download the corresponding paired end reads from SRA. If both local files and SRA accession are available, the local files will be preferred.

Missing values can be specified by empty columns or by writing `NA`.

##### Calling scenario

Varlociraptor supports integrated uncertainty aware calling and filtering of variants for arbitrary scenarios. These are defined as so-called scenarios, via a [variant calling grammar](https://varlociraptor.github.io/docs/calling#generic-variant-calling).
* For each group, a scenario is rendered via [Jinja](https://jinja.palletsprojects.com).
* Therefore, edit the template scenario (`scenario.yaml`) according to your needs. The sample sheet is available for jinja rendering as a pandas data frame in the variable `samples`. This allows to customize the scenario according to the contents of the sample sheet. You can therefore add additional columns to the sample sheet (e.g. purity) and access them in the scenario template, in order to pass the information to Varlociraptor.

#### Step 3: Execute workflow

Test your configuration by performing a dry-run via

    snakemake --use-conda -n

Execute the workflow locally via

    snakemake --use-conda --cores $N --set-scatter calls=$M

using `$N` cores and `$M` separate calling and annotation jobs per sample.
Alternatively, the workflow can be e.g. executed in a cluster environment via

    snakemake --use-conda --cluster qsub --jobs 100 --set-scatter calls=$M

or

    snakemake --use-conda --drmaa --jobs 100 --set-scatter calls=$M

If you not only want to fix the software stack but also the underlying OS, use

    snakemake --use-conda --use-singularity

in combination with any of the modes above.

Naturally, Snakemake workflow can also be executed **in the cloud**, and there are **specialized profiles available for all major cluster engines**.
See the [Snakemake documentation](https://snakemake.readthedocs.io/en/stable/executable.html) for further details.

#### Step 4: Investigate results

After successful execution, you can create a self-contained interactive HTML report with all results via:

    snakemake --report report.zip

This report can, e.g., be forwarded to your collaborators.

#### Step 5: Commit changes

Whenever you change something, don't forget to commit the changes back to your github copy of the repository:

    git commit -a
    git push

#### Step 6: Obtain updates from upstream

Whenever you want to synchronize your workflow copy with new developments from upstream, do the following.

1. Once, register the upstream repository in your local copy: `git remote add -f upstream git@github.com:snakemake-workflows/dna-seq-varlociraptor.git` or `git remote add -f upstream https://github.com/snakemake-workflows/dna-seq-varlociraptor.git` if you do not have setup ssh keys.
2. Update the upstream version: `git fetch upstream`.
3. Create a diff with the current version: `git diff HEAD upstream/master workflow .test > upstream-changes.diff`.
4. Investigate the changes: `vim upstream-changes.diff`.
5. Apply the modified diff via: `git apply upstream-changes.diff`.
6. Carefully check whether you need to update the config files: `git diff HEAD upstream/master config`. If so, do it manually, and only where necessary, since you would otherwise likely overwrite your settings and samples.


#### Step 7: Contribute back

In case you have also changed or added steps, please consider contributing them back to the original repository:

1. [Fork](https://help.github.com/en/articles/fork-a-repo) the original repo to a personal or lab account.
2. [Clone](https://help.github.com/en/articles/cloning-a-repository) the fork to your local system, to a different place than where you ran your analysis.
3. Copy the modified files from your analysis to the clone of your fork, e.g., `cp -r workflow path/to/fork`. Make sure to **not** accidentally copy config file contents or sample sheets. Instead, manually update the example config files if necessary.
4. Commit and push your changes to your fork.
5. Create a [pull request](https://help.github.com/en/articles/creating-a-pull-request) against the original repository.

## Testing

Test cases are in the subfolder `.test`. They are automtically executed via continuous integration with Github actions.
