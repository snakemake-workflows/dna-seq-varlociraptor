from snakemake.utils import min_version

min_version("5.7.0")

import pandas as pd
from snakemake.remote import FTP

ftp = FTP.RemoteProvider()

samples = pd.read_csv(config["samples"], sep="\t").set_index("sample_name", drop=False).sort_index()

def _group_or_sample(row):
    group = row.get("group", None)
    if pd.isnull(group):
        return row["sample_name"]
    return group

samples["group"] = [_group_or_sample(row) for _, row in samples.iterrows()]

units = pd.read_csv(config["units"], sep="\t", dtype=str).set_index(["sample_name", "unit_name"], drop=False).sort_index()

def is_paired_end(sample):
    assert not units.loc[sample, "fq1"].isnull().any()
    fq2_null = units.loc[sample, "fq2"].isnull()
    all_paired = not fq2_null.any()
    assert fq2_null.all() or all_paired, "invalid units for sample {}, must be all paired end or all single end".format(sample)
    return all_paired


def get_merged(wildcards):
    if is_paired_end(wildcards.sample):
        return ["merged/{sample}.1.fastq.gz",
                "merged/{sample}.2.fastq.gz"]
    return "merged/{sample}.fastq.gz"

def get_group_aliases(wildcards):
    return samples.loc[samples["group"] == wildcards.group]["alias"]


def get_group_samples(wildcards):
    return samples.loc[samples["group"] == wildcards.group]["sample_name"]


def get_group_bams(wildcards):
    return expand("recal/{sample}.sorted.bam", sample=get_group_samples(wildcards))

def get_group_observations(wildcards):
    return expand("observations/{group}/{sample}.{caller}.bcf", 
                  caller=wildcards.caller, 
                  group=wildcards.group,
                  sample=get_group_samples(wildcards))

def get_group_bais(wildcards):
    return expand("recal/{sample}.sorted.bam.bai", sample=get_group_samples(wildcards))

def is_activated(xpath):
    c = config
    for entry in xpath.split("/"):
        c = c.get(entry, {})
    return bool(c["activate"])


def get_read_group(wildcards):
    """Denote sample name and platform in read group."""
    return r"-R '@RG\tID:{sample}\tSM:{sample}\tPL:{platform}'".format(
        sample=wildcards.sample,
        platform=samples.loc[wildcards.sample, "platform"])


def get_tmb_targets():
    if is_activated("tmb"):
        return expand("plots/tmb/{group}.tmb.svg",
                      group=groups)
    else:
        return []


def get_annotated_bcf(wildcards, group=None):
    if group is None:
        group = wildcards.group
    selection = ".annotated"
    if is_activated("annotations/vcfs"):
        selection += ".db-annotated"
    if is_activated("annotations/dbnsfp"):
        selection += ".dbnsfp"
    if is_activated("annotations/dgidb"):
        selection += ".dgidb"
    return "calls/{group}{selection}.bcf".format(group=group, selection=selection)


wildcard_constraints:
    group="|".join(samples["group"].unique()),
    sample="|".join(samples["sample_name"]),
    caller="|".join(["freebayes", "delly"])

caller=list(filter(None, ["freebayes" if is_activated("calling/freebayes") else None, "delly" if is_activated("calling/delly") else None]))
