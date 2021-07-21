import glob
from os import path

import pandas as pd
from snakemake.remote import FTP
from snakemake.utils import validate

ftp = FTP.RemoteProvider()

validate(config, schema="../schemas/config.schema.yaml")

samples = (
    pd.read_csv(
        config["samples"],
        sep="\t",
        dtype={"sample_name": str, "group": str},
        comment="#",
    )
    .set_index("sample_name", drop=False)
    .sort_index()
)


def get_final_output():
    final_output = []

    if config["report"]["activate"]:
        final_output.extend(
            expand(
                "results/vcf-report/all.{event}/",
                event=config["calling"]["fdr-control"]["events"],
            )
        )
    else:
        final_output.extend(
            expand(
                "results/final-calls/{group}.{event}.fdr-controlled.bcf",
                group=groups,
                event=config["calling"]["fdr-control"]["events"],
            )
        )

    if config["tables"]["activate"]:
        final_output.extend(
            expand(
                "results/tables/{group}.{event}.fdr-controlled.tsv",
                group=groups,
                event=config["calling"]["fdr-control"]["events"],
            )
        )
        if config["tables"].get("generate_excel", False):
            final_output.extend(
                expand(
                    "results/tables/{group}.{event}.fdr-controlled.xlsx",
                    group=groups,
                    event=config["calling"]["fdr-control"]["events"],
                )
            )

    final_output.extend(get_mutational_burden_targets())

    return final_output


def _group_or_sample(row):
    group = row.get("group", None)
    if pd.isnull(group):
        return row["sample_name"]
    return group


samples["group"] = [_group_or_sample(row) for _, row in samples.iterrows()]
validate(samples, schema="../schemas/samples.schema.yaml")

units = (
    pd.read_csv(
        config["units"],
        sep="\t",
        dtype={"sample_name": str, "unit_name": str},
        comment="#",
    )
    .set_index(["sample_name", "unit_name"], drop=False)
    .sort_index()
)
validate(units, schema="../schemas/units.schema.yaml")


def get_gather_calls_input(ext="bcf"):
    def inner(wildcards):
        if wildcards.by == "odds":
            pattern = "results/calls/{{{{group}}}}.{{{{event}}}}.{{{{filter}}}}.{{scatteritem}}.filtered_odds.{ext}"
        elif wildcards.by == "ann":
            pattern = "results/calls/{{{{group}}}}.{{{{filter}}}}.{{scatteritem}}.filtered_ann.{ext}"
        else:
            raise ValueError(
                "Unexpected wildcard value for 'by': {}".format(wildcards.by)
            )
        return gather.calling(pattern.format(ext=ext))

    return inner


def get_control_fdr_input(wildcards):
    query = get_fdr_control_params(wildcards)
    if not is_activated("benchmarking"):
        by = "ann" if query["local"] else "odds"
        return "results/calls/{{group}}.{{event}}.{{filter}}.filtered_{by}.bcf".format(
            by=by
        )
    else:
        return "results/calls/{group}.bcf"


def get_recalibrate_quality_input(wildcards, bai=False):
    ext = "bai" if bai else "bam"
    if is_activated("calc_consensus_reads"):
        return "results/consensus/{}.sorted.{}".format(wildcards.sample, ext)
    elif is_activated("remove_duplicates"):
        return "results/dedup/{}.sorted.{}".format(wildcards.sample, ext)
    else:
        return "results/mapped/{}.sorted.{}".format(wildcards.sample, ext)


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
        return "pipe/cutadapt/{S}/{U}.fq1.fastq{E}".format(
            S=unit.sample_name, U=unit.unit_name, E=ending
        )
    else:
        # paired end local sample
        return expand(
            "pipe/cutadapt/{S}/{U}.{{read}}.fastq{E}".format(
                S=unit.sample_name, U=unit.unit_name, E=ending
            ),
            read=["fq1", "fq2"],
        )


def get_cutadapt_pipe_input(wildcards):
    pattern = units.loc[wildcards.sample].loc[wildcards.unit, wildcards.fq]
    if "*" in pattern:
        files = sorted(
            glob.glob(units.loc[wildcards.sample].loc[wildcards.unit, wildcards.fq])
        )
        if not files:
            raise ValueError(
                "No raw fastq files found for unit pattern {} (sample {}). "
                "Please check the your sample sheet.".format(
                    wildcards.unit, wildcards.sample
                )
            )
    else:
        files = [pattern]

    return files


def get_cutadapt_adapters(wildcards):
    unit = units.loc[wildcards.sample].loc[wildcards.unit]
    try:
        adapters = unit["adapters"]
        if isinstance(adapters, str):
            return adapters
        return ""
    except KeyError:
        return ""


def is_paired_end(sample):
    sample_units = units.loc[sample]
    fq2_null = sample_units["fq2"].isnull()
    sra_null = sample_units["sra"].isnull()
    paired = ~fq2_null | ~sra_null
    all_paired = paired.all()
    all_single = (~paired).all()
    assert (
        all_single or all_paired
    ), "invalid units for sample {}, must be all paired end or all single end".format(
        sample
    )
    return all_paired


def group_is_paired_end(group):
    samples = get_group_samples(group)
    return all([is_paired_end(sample) for sample in samples])


def get_map_reads_input(wildcards):
    if is_paired_end(wildcards.sample):
        return [
            "results/merged/{sample}_R1.fastq.gz",
            "results/merged/{sample}_R2.fastq.gz",
        ]
    return "results/merged/{sample}_single.fastq.gz"


def get_group_aliases(wildcards):
    return samples.loc[samples["group"] == wildcards.group]["alias"]


def get_group_samples(group):
    return samples.loc[samples["group"] == group]["sample_name"]


def get_group_sample_aliases(wildcards, controls=True):
    if controls:
        return samples.loc[samples["group"] == wildcards.group]["alias"]
    return samples.loc[
        (samples["group"] == wildcards.group) & (samples["control"] == "no")
    ]["alias"]


def get_sample_bam(wildcards, bai=False):
    ext = "bai" if bai else "bam"
    if is_activated("primers/trimming"):
        if group_is_paired_end(wildcards.group):
            return "results/trimmed/{sample}.trimmed.{ext}".format(
                sample=wildcards.sample, ext=ext
            )
        else:
            WorkflowError("Primer trimming is only available for paired end data.")
    return "results/recal/{sample}.sorted.{ext}".format(
        sample=wildcards.sample, ext=ext
    )


def get_group_bams(wildcards, bai=False):
    ext = "bai" if bai else "bam"
    if is_activated("primers/trimming"):
        if group_is_paired_end(wildcards.group):
            return expand(
                "results/trimmed/{sample}.trimmed.{ext}",
                sample=get_group_samples(wildcards.group),
                ext=ext,
            )
        else:
            WorkflowError("Primer trimming is only available for paired end data.")
    return expand(
        "results/recal/{sample}.sorted.{ext}",
        sample=get_group_samples(wildcards.group),
        ext=ext,
    )


def get_group_bams_report(group):
    if is_activated("primers/trimming"):
        if group_is_paired_end(group):
            return [
                (sample, "results/trimmed/{}.trimmed.bam".format(sample))
                for sample in get_group_samples(group)
            ]
        else:
            WorkflowError("Primer trimming is only available for paired end data.")
    return [
        (sample, "results/recal/{}.sorted.bam".format(sample))
        for sample in get_group_samples(group)
    ]


def _get_batch_info(wildcards, yield_sample=False, yield_bam=False):
    for group in get_report_batch(wildcards):
        for sample, bam in get_group_bams_report(group):
            if yield_sample and yield_bam:
                yield sample, bam
            elif yield_sample:
                yield sample
            elif yield_bam:
                yield bam
            else:
                raise ValueError("Either set yield_sample or yield_bam.")


def get_batch_bams(wildcards):
    yield from _get_batch_info(wildcards, yield_bam=True)


def get_report_bam_params(wildcards, input):
    return [
        "{group}:{sample}={bam}".format(group=wildcards.group, sample=sample, bam=bam)
        for sample, bam in zip(
            _get_batch_info(wildcards, yield_sample=True), input.bams
        )
    ]


def get_batch_bcfs(wildcards):
    for group in get_report_batch(wildcards):
        yield "results/final-calls/{group}.{event}.fdr-controlled.bcf".format(
            group=group, event=wildcards.event
        )


def get_report_bcf_params(wildcards, input):
    return [
        "{group}={bcf}".format(group=group, bcf=bcf)
        for group, bcf in zip(get_report_batch(wildcards), input.bcfs)
    ]


def get_consensus_input(wildcards):
    if wildcards.read_type == "se":
        return "results/consensus/fastq/{}.se.fq".format(wildcards.sample)
    return [
        "results/consensus/fastq/{}.1.fq".format(wildcards.sample),
        "results/consensus/fastq/{}.2.fq".format(wildcards.sample),
    ]


def get_resource(name):
    return workflow.source_path("../resources/{}".format(name))


# TODO Can be reduced to "results/regions/{group}.target_regions.bed" when single primer trimming gets implemented
def get_regions(wildcards):
    if is_activated("primers/trimming"):
        return "results/primers/target_regions.merged.bed"
    return "results/regions/{group}.target_regions.filtered.bed".format(
        group=wildcards.group
    )


# TODO Can be reduced to "results/regions/{group}.excluded_regions.bed" when single primer trimming gets implemented
def get_excluded_regions(wildcards):
    if is_activated("primers/trimming"):
        return "results/primers/excluded_regions.bed"
    return "results/regions/{group}.excluded_regions.bed".format(group=wildcards.group)


def get_group_observations(wildcards):
    # TODO if group contains only a single sample, do not require sorting.
    return expand(
        "results/observations/{group}/{sample}.{caller}.{scatteritem}.bcf",
        caller=wildcards.caller,
        group=wildcards.group,
        scatteritem=wildcards.scatteritem,
        sample=get_group_samples(wildcards.group),
    )


def is_activated(xpath):
    c = config
    for entry in xpath.split("/"):
        c = c.get(entry, {})
    return bool(c.get("activate", False))


def get_read_group(wildcards):
    """Denote sample name and platform in read group."""
    return r"-R '@RG\tID:{sample}\tSM:{sample}\tPL:{platform}'".format(
        sample=wildcards.sample, platform=samples.loc[wildcards.sample, "platform"]
    )


def get_mutational_burden_targets():
    if is_activated("mutational_burden"):
        return expand(
            "results/plots/mutational-burden/{sample.group}.{sample.sample_name}.{mode}.mutational-burden.svg",
            mode=config["mutational_burden"].get("mode", "curve"),
            sample=samples.itertuples(),
        )
    else:
        return []


def get_mutational_burden_events(wildcards):
    try:
        events = samples.loc[wildcards.sample, "mutational_burden_events"]
    except KeyError:
        events = None
    if pd.isna(events):
        events = config["mutational_burden"]["events"]
    else:
        events = map(str.strip, events.split(","))
    return " ".join(events)


def get_scattered_calls(ext="bcf"):
    def inner(wildcards):
        return expand(
            "results/calls/{{group}}.{caller}.{{scatteritem}}.sorted.{ext}",
            caller=caller,
            ext=ext,
        )

    return inner


def get_selected_annotations():
    selection = ".annotated"
    if is_activated("annotations/vcfs"):
        selection += ".db-annotated"
    if is_activated("annotations/dgidb"):
        selection += ".dgidb"
    return selection


def get_annotated_bcf(wildcards):
    selection = get_selected_annotations()
    return "results/calls/{group}.{scatteritem}{selection}.bcf".format(
        group=wildcards.group, selection=selection, scatteritem=wildcards.scatteritem
    )


def get_gather_annotated_calls_input(ext="bcf"):
    def inner(wildcards):
        selection = get_selected_annotations()
        return gather.calling(
            "results/calls/{{{{group}}}}.{{scatteritem}}{selection}.{ext}".format(
                ext=ext, selection=selection
            )
        )

    return inner


def get_candidate_calls():
    filter = config["calling"]["filter"].get("candidates")
    if filter:
        return "results/candidate-calls/{group}.{caller}.{scatteritem}.filtered.bcf"
    else:
        return "results/candidate-calls/{group}.{caller}.{scatteritem}.bcf"


def get_report_batch(wildcards):
    if wildcards.batch == "all":
        groups = samples["group"].unique()
    else:
        groups = samples.loc[
            samples[config["report"]["stratify"]["by-column"]] == wildcards.batch,
            "group",
        ].unique()
    if not any(groups):
        raise ValueError("No samples found. Is your sample sheet empty?")
    return groups


def get_merge_calls_input(ext="bcf"):
    def inner(wildcards):
        return expand(
            "results/calls/{{group}}.{vartype}.{{event}}.{filter}.fdr-controlled.{ext}",
            ext=ext,
            vartype=["SNV", "INS", "DEL", "MNV", "BND", "INV", "DUP", "REP"],
            filter=config["calling"]["fdr-control"]["events"][wildcards.event][
                "filter"
            ],
        )

    return inner


def get_merge_calls_input_report(wildcards, ext="bcf"):
    return expand(
        "{{group}}.{vartype}.{{event}}.{filter}=results/calls/{{group}}.{vartype}.{{event}}.{filter}.fdr-controlled.{ext}",
        ext=ext,
        vartype=["SNV", "INS", "DEL", "MNV", "BND", "INV", "DUP", "REP"],
        filter=config["calling"]["fdr-control"]["events"][wildcards.event]["filter"],
    )


def get_vep_threads():
    n = len(samples)
    if n:
        return max(workflow.cores / n, 1)
    else:
        return 1


def get_fdr_control_params(wildcards):
    query = config["calling"]["fdr-control"]["events"][wildcards.event]
    threshold = query.get(
        "threshold", config["calling"]["fdr-control"].get("threshold", 0.05)
    )
    events = query["varlociraptor"]
    local = (
        "--local"
        if query.get("local", config["calling"]["fdr-control"].get("local", False))
        else ""
    )
    return {"threshold": threshold, "events": events, "local": local}


def get_fixed_candidate_calls(wildcards):
    if wildcards.caller == "delly":
        return "results/candidate-calls/{group}.delly.no_bnds.bcf"
    else:
        return "results/candidate-calls/{group}.{caller}.bcf"


wildcard_constraints:
    group="|".join(samples["group"].unique()),
    sample="|".join(samples["sample_name"]),
    caller="|".join(["freebayes", "delly"]),
    filter="|".join(config["calling"]["filter"]),


caller = list(
    filter(
        None,
        [
            "freebayes" if is_activated("calling/freebayes") else None,
            "delly" if is_activated("calling/delly") else None,
        ],
    )
)

###### Annotations ########

annotations = [
    (e, f) for e, f in config["annotations"]["vcfs"].items() if e != "activate"
]


def get_annotation_pipes(wildcards, input):
    if annotations:
        return "| {}".format(
            " | ".join(
                [
                    "SnpSift annotate -name {prefix}_ {path} /dev/stdin".format(
                        prefix=prefix, path=path
                    )
                    for (prefix, _), path in zip(annotations, input.annotations)
                ]
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
    return expand(
        "results/trimmed/{sample}/{unit}_{read}.fastq.gz",
        unit=units.loc[wc.sample, "unit_name"],
        sample=wc.sample,
        read=wc.read,
    )


def get_vembrane_expression(wc):
    expression = (
        config["tables"]
        .get("output", {})
        .get(
            "expression",
            "INDEX, CHROM, POS, REF, ALT[0], ANN['Consequence'], ANN['IMPACT'], ANN['SYMBOL'], ANN['Feature']",
        )
    )
    parts = [expression]
    if config["tables"].get("output", {}).get("event_prob", False):
        parts.append(
            ", ".join(
                f"1-10**(-INFO['PROB_{x.upper()}']/10)"
                for x in config["calling"]["fdr-control"]["events"][wc.event][
                    "varlociraptor"
                ]
            )
        )
    if config["tables"].get("output", {}).get("genotype", False):
        parts.append(
            ", ".join(
                f"FORMAT['AF']['{sample}']" for sample in get_group_sample_aliases(wc)
            )
        )
    if config["tables"].get("output", {}).get("depth", False):
        parts.append(
            ", ".join(
                f"FORMAT['DP']['{sample}']" for sample in get_group_sample_aliases(wc)
            )
        )
    return ", ".join(parts)
