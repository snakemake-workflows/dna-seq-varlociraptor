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

primer_panels = (
    (
        pd.read_csv(
            config["primers"]["trimming"]["tsv"],
            sep="\t",
            dtype={"panel": str, "fa1": str, "fa2": str},
            comment="#",
        )
        .set_index(["panel"], drop=False)
        .sort_index()
    )
    if config["primers"]["trimming"].get("tsv", "")
    else None
)


def get_gather_calls_input(ext="bcf"):
    def inner(wildcards):
        if wildcards.by == "odds":
            pattern = "results/calls/{{{{group}}}}.{{{{event}}}}.{{scatteritem}}.filtered_odds.{ext}"
        elif wildcards.by == "ann":
            pattern = "results/calls/{{{{group}}}}.{{{{event}}}}.{{scatteritem}}.filtered_ann.{ext}"
        else:
            raise ValueError(
                "Unexpected wildcard value for 'by': {}".format(wildcards.by)
            )
        return gather.calling(pattern.format(ext=ext))

    return inner


def get_control_fdr_input(wildcards):
    query = get_fdr_control_params(wildcards)
    if not is_activated("benchmarking") and query["filter"]:
        by = "ann" if query["local"] else "odds"
        return "results/calls/{{group}}.{{event}}.filtered_{by}.bcf".format(by=by)
    else:
        return "results/calls/{group}.bcf"


def get_recalibrate_quality_input(wildcards, bai=False):
    ext = "bai" if bai else "bam"
    if is_activated("calc_consensus_reads"):
        return "results/consensus/{}.sorted.{}".format(wildcards.sample, ext)
    elif is_activated("primers/trimming"):
        return "results/trimmed/{sample}.trimmed.{ext}".format(
            sample=wildcards.sample, ext=ext
        )
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


def get_consensus_input(wildcards):
    if is_activated("primers/trimming"):
        return "results/trimmed/{}.trimmed.bam".format(wildcards.sample)
    elif is_activated("remove_duplicates"):
        return "results/dedup/{}.sorted.bam".format(wildcards.sample)
    else:
        return "results/mapped/{}.sorted.bam".format(wildcards.sample)


def get_trimming_input(wildcards):
    if is_activated("remove_duplicates"):
        return "results/dedup/{}.sorted.bam".format(wildcards.sample)
    else:
        return "results/mapped/{}.sorted.bam".format(wildcards.sample)


def get_primer_bed(wc):
    if isinstance(primer_panels, pd.DataFrame):
        if not pd.isna(primer_panels.loc[wc.panel, "fa2"]):
            return "results/primers/{}_primers.bedpe".format(wc.panel)
        else:
            return "results/primers/{}_primers.bed".format(wc.panel)
    else:
        if config["primers"]["trimming"].get("primers_fa2", ""):
            return "results/primers/uniform_primers.bedpe"
        else:
            return "results/primers/uniform_primers.bed"


def get_sample_primer_fastas(sample):
    if isinstance(primer_panels, pd.DataFrame):
        panel = samples.loc[sample, "panel"]
        if not pd.isna(primer_panels.loc[panel, "fa2"]):
            return [
                primer_panels.loc[panel, "fa1"],
                primer_panels.loc[panel, "fa2"],
            ]
        return primer_panels.loc[panel, "fa1"]
    else:
        if config["primers"]["trimming"].get("primers_fa2", ""):
            return [
                config["primers"]["trimming"]["primers_fa1"],
                config["primers"]["trimming"]["primers_fa2"],
            ]
        return config["primers"]["trimming"]["primers_fa1"]


def get_panel_primer_input(panel):
    if panel == "uniform":
        if config["primers"]["trimming"].get("primers_fa2", ""):
            return [
                config["primers"]["trimming"]["primers_fa1"],
                config["primers"]["trimming"]["primers_fa2"],
            ]
        return config["primers"]["trimming"]["primers_fa1"]
    else:
        panel = primer_panels.loc[panel]
        if not pd.isna(panel["fa2"]):
            return [panel["fa1"], panel["fa2"]]
        return panel["fa1"]


def input_is_fasta(primers):
    primers = primers[0] if isinstance(primers, list) else primers
    fasta_suffixes = ("fasta", "fa")
    return True if primers.endswith(fasta_suffixes) else False


def get_primer_regions(wc):
    if isinstance(primer_panels, pd.DataFrame):
        return "results/primers/{}_primer_regions.tsv".format(
            samples.loc[wc.sample, "panel"]
        )
    return "results/primers/uniform_primer_regions.tsv"


def get_markduplicates_extra(wc):
    c = config["params"]["picard"]["MarkDuplicates"]

    if units.loc[wc.sample]["umis"].isnull().any():
        b = ""
    else:
        b = "--BARCODE_TAG RX"

    if is_activated("calc_consensus_reads"):
        d = "--TAG_DUPLICATE_SET_MEMBERS true"
    else:
        d = ""

    return f"{c} {b} {d}"


def get_group_bams(wildcards, bai=False):
    ext = "bai" if bai else "bam"
    if is_activated("primers/trimming") and not group_is_paired_end(wildcards.group):
        WorkflowError("Primer trimming is only available for paired end data.")
    return expand(
        "results/recal/{sample}.sorted.{ext}",
        sample=get_group_samples(wildcards.group),
        ext=ext,
    )


def get_group_bams_report(group):
    return [
        (sample, "results/recal/{}.sorted.bam".format(sample))
        for sample in get_group_samples(group)
    ]


def _get_batch_info(wildcards):
    for group in get_report_batch(wildcards):
        for sample, bam in get_group_bams_report(group):
            yield sample, bam, group


def get_batch_bams(wildcards):
    return (bam for _, bam, _ in _get_batch_info(wildcards))


def get_report_bam_params(wildcards, input):
    return [
        "{group}:{sample}={bam}".format(group=group, sample=sample, bam=bam)
        for (sample, _, group), bam in zip(_get_batch_info(wildcards), input.bams)
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


def get_processed_consensus_input(wildcards):
    if wildcards.read_type == "se":
        return "results/consensus/fastq/{}.se.fq".format(wildcards.sample)
    return [
        "results/consensus/fastq/{}.1.fq".format(wildcards.sample),
        "results/consensus/fastq/{}.2.fq".format(wildcards.sample),
    ]


def get_resource(name):
    return workflow.source_path("../resources/{}".format(name))


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
            "results/calls/{{group}}.{vartype}.{{event}}.fdr-controlled.{ext}",
            ext=ext,
            vartype=["SNV", "INS", "DEL", "MNV", "BND", "INV", "DUP", "REP"],
        )

    return inner


def get_vep_threads():
    n = len(samples)
    if n:
        return max(workflow.cores / n, 1)
    else:
        return 1


def get_plugin_aux(plugin, index=False):
    if plugin in config["annotations"]["vep"]["plugins"]:
        if plugin == "REVEL":
            suffix = ".tbi" if index else ""
            return "resources/revel_scores.tsv.gz{suffix}".format(suffix=suffix)
    return []


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
    return {
        "threshold": threshold,
        "events": events,
        "local": local,
        "filter": query.get("filter"),
    }


def get_fixed_candidate_calls(wildcards):
    if wildcards.caller == "delly":
        return "results/candidate-calls/{group}.delly.no_bnds.bcf"
    else:
        return "results/candidate-calls/{group}.{caller}.bcf"


def get_filter_targets(wildcards, input):
    if input.predefined:
        return " | bedtools intersect -a /dev/stdin -b {input.predefined} ".format(
            input=input
        )
    else:
        return ""


def get_annotation_filter(wildcards):
    filter = config["calling"]["fdr-control"]["events"][wildcards.event]["filter"]
    filter = (
        [config["calling"]["filter"][filter]]
        if isinstance(filter, str)
        else map(lambda x: config["calling"]["filter"][x], filter)
    )
    return " and ".join(filter)


wildcard_constraints:
    group="|".join(samples["group"].unique()),
    sample="|".join(samples["sample_name"]),
    caller="|".join(["freebayes", "delly"]),
    filter="|".join(config["calling"]["filter"]),
    event="|".join(config["calling"]["fdr-control"]["events"].keys()),


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


def get_tabix_revel_params():
    # Indexing of REVEL-score file where the column depends on the reference
    column = 2 if config["ref"]["build"] == "GRCh37" else 3
    return f"-f -s 1 -b {column} -e {column}"


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
                f"10**(-INFO['PROB_{x.upper()}']/10)"
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


def get_sample_alias(wildcards):
    return samples.loc[wildcards.sample, "alias"]


def get_dgidb_datasources():
    if config["annotations"]["dgidb"].get("datasources", ""):
        return "-s {}".format(" ".join(config["annotations"]["dgidb"]["datasources"]))
    return ""


def get_bowtie_insertsize():
    if config["primers"]["trimming"].get("library_length", 0) != 0:
        return "-X {}".format(config["primers"]["trimming"].get("library_length"))
    return ""


def get_filter_params(wc):
    if isinstance(get_panel_primer_input(wc.panel), list):
        return "-b -f 2"
    return "-b -F 4"


def get_single_primer_flag(wc):
    if not isinstance(get_sample_primer_fastas(wc.sample), list):
        return "--first-of-pair"
    return ""


def format_bowtie_primers(wc, primers):
    if isinstance(primers, list):
        return "-1 {r1} -2 {r2}".format(r1=primers[0], r2=primers[1])
    return primers
