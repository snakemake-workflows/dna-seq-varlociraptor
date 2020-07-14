import glob

import pandas as pd
from snakemake.remote import FTP
from snakemake.utils import validate

ftp = FTP.RemoteProvider()

validate(config, schema="../schemas/config.schema.yaml")

samples = pd.read_csv(config["samples"], sep="\t", dtype={"sample_name": str, "group": str}).set_index("sample_name", drop=False).sort_index()

def get_final_output():
    if config["igv_report"]["activate"]:
        final_output = expand("results/igv-report/{group}.{event}.html",
                        group=groups,
                        event=config["calling"]["fdr-control"]["events"]),
    else:
        final_output = expand("results/merged-calls/{group}.{event}.fdr-controlled.bcf",
                        group=groups,
                        event=config["calling"]["fdr-control"]["events"]),
    return final_output

def _group_or_sample(row):
    group = row.get("group", None)
    if pd.isnull(group):
        return row["sample_name"]
    return group

samples["group"] = [_group_or_sample(row) for _, row in samples.iterrows()]
validate(samples, schema="../schemas/samples.schema.yaml")

units = pd.read_csv(config["units"], sep="\t", dtype={"sample_name": str, "unit_name": str}).set_index(["sample_name", "unit_name"], drop=False).sort_index()
validate(units, schema="../schemas/units.schema.yaml")


def get_recalibrate_quality_input(wildcards, bai=False):
    ext = ".bai" if bai else ""
    if is_activated("remove_duplicates"):
        return "results/dedup/{}.sorted.bam{}".format(wildcards.sample, ext)
    else:
        return "results/mapped/{}.sorted.bam{}".format(wildcards.sample, ext)


def get_cutadapt_input(wildcards):
    unit = units.loc[wildcards.sample].loc[wildcards.unit]

    if pd.isna(unit["fq1"]):
        # SRA sample (always paired-end for now)
        accession = unit["sra"]
        return expand("sra/{accession}_{read}.fastq", accession=accession, read=[1, 2])

    if unit["fq1"].endswith("gz"):
        ending = ".gz"
    else:
        ending = ""

    if pd.isna(unit["fq2"]):
        # single end local sample
        return "pipe/cutadapt/{S}/{U}.fq1.fastq{E}".format(S=unit.sample_name, U=unit.unit_name, E=ending)
    else:
        # paired end local sample
        return expand("pipe/cutadapt/{S}/{U}.{{read}}.fastq{E}".format(S=unit.sample_name, U=unit.unit_name, E=ending), read=["fq1","fq2"])


def get_cutadapt_pipe_input(wildcards):
    files = list(sorted(glob.glob(units.loc[wildcards.sample].loc[wildcards.unit, wildcards.fq])))
    assert(len(files) > 0)
    return files


def is_paired_end(sample):
    sample_units = units.loc[sample]
    fq2_null = sample_units["fq2"].isnull()
    sra_null = sample_units["sra"].isnull()
    paired = ~fq2_null | ~sra_null
    all_paired = paired.all()
    all_single = (~paired).all()
    assert all_single or all_paired, "invalid units for sample {}, must be all paired end or all single end".format(sample)
    return all_paired

def group_is_paired_end(wildcards):
    samples = get_group_samples(wildcards)
    return all([is_paired_end(sample) for sample in samples])

def get_map_reads_input(wildcards):
    if is_paired_end(wildcards.sample):
        return ["results/merged/{sample}_R1.fastq.gz",
                "results/merged/{sample}_R2.fastq.gz"]
    return "results/merged/{sample}_single.fastq.gz"


def get_group_aliases(wildcards):
    return samples.loc[samples["group"] == wildcards.group]["alias"]


def get_group_samples(wildcards):
    return samples.loc[samples["group"] == wildcards.group]["sample_name"]


def get_group_bams(wildcards, bai=False):
    ext = ".bai" if bai else ""
    if group_is_paired_end(wildcards) and is_activated("primers/trimming"):
        return expand("results/trimmed/{sample}.trimmed.bam{ext}", sample=get_group_samples(wildcards), ext=ext)
    return expand("results/recal/{sample}.sorted.bam{ext}", sample=get_group_samples(wildcards), ext=ext)


def get_regions():
    if is_activated("primers/trimming"):
        return "results/primers/target_regions.merged.bed"
    else:
        return []

def get_excluded_regions():
    if is_activated("primers/trimming"):
        return "results/primers/excluded_regions.bed"
    else:
        return []

def get_group_observations(wildcards):
    return expand("results/observations/{group}/{sample}.{caller}.bcf", 
                  caller=wildcards.caller, 
                  group=wildcards.group,
                  sample=get_group_samples(wildcards))

def is_activated(xpath):
    c = config
    for entry in xpath.split("/"):
        c = c.get(entry, {})
    return bool(c.get("activate", False))

def get_read_group(wildcards):
    """Denote sample name and platform in read group."""
    return r"-R '@RG\tID:{sample}\tSM:{sample}\tPL:{platform}'".format(
        sample=wildcards.sample,
        platform=samples.loc[wildcards.sample, "platform"])


def get_tmb_targets():
    if is_activated("tmb"):
        return expand("results/plots/tmb/{group}.{mode}.tmb.svg",
                      group=groups,
                      mode=config["tmb"].get("mode", "curve"))
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
    return "results/calls/{group}{selection}.bcf".format(group=group, selection=selection)


def get_oncoprint_batch(wildcards):
    if wildcards.batch == "all":
        groups = samples["group"].unique()
    else:
        groups = samples.loc[samples[config["oncoprint"]["stratify"]["by-column"]] == wildcards.batch, "group"].unique()
    if not any(groups):
        raise ValueError("No samples found. Is your sample sheet empty?")
    return groups


def get_merge_calls_input(ext=".bcf"):
    def inner(wildcards):
        return expand("results/calls/{{group}}.{vartype}.{{event}}.{filter}.fdr-controlled{ext}",
                      ext=ext,
                      vartype=["SNV", "INS", "DEL", "MNV"],
                      filter=config["calling"]["fdr-control"]["events"][wildcards.event]["filter"])
    return inner


wildcard_constraints:
    group="|".join(samples["group"].unique()),
    sample="|".join(samples["sample_name"]),
    caller="|".join(["freebayes", "delly"])

caller=list(filter(None, ["freebayes" if is_activated("calling/freebayes") else None, "delly" if is_activated("calling/delly") else None]))

###### Annotations ########

annotations = [(e, f) for e, f in config["annotations"]["vcfs"].items() if e != "activate"]

def get_annotation_pipes(wildcards, input):
     if annotations:
         return "| {}".format(" | ".join(
             ["SnpSift annotate -name {prefix}_ {path} /dev/stdin".format(prefix=prefix, path=path)
              for (prefix, _), path in zip(annotations, input.annotations)]
              )
         )
     else:
         return ""


def get_annotation_vcfs(idx=False):
    fmt = lambda f: f if not idx else "{}.tbi".format(f)
    return [fmt(f) for _, f in annotations]


def get_tabix_params(wildcards):
    if wildcards.format == "vcf":
        return "-p vcf"
    if wildcards.format == "txt":
        return "-s 1 -b 2 -e 2"
    raise ValueError("Invalid format for tabix: {}".format(wildcards.format))


def get_fastqs(wc):
    return expand("results/trimmed/{sample}/{unit}_{read}.fastq.gz", unit=units.loc[wc.sample, "unit_name"], sample=wc.sample, read=wc.read)
