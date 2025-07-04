import glob
from os import path

import yaml
import pandas as pd
from snakemake.utils import validate

validate(config, schema="../schemas/config.schema.yaml")

samples = (
    pd.read_csv(
        config["samples"],
        sep="\t",
        dtype={"sample_name": str, "group": str, "umi_len": "Int64"},
        comment="#",
    )
    .set_index("sample_name", drop=False)
    .sort_index()
)
if not "mutational_burden_events" in samples.columns:
    samples["mutational_burden_events"] = pd.NA

# construct genome name
datatype_genome = "dna"
species = config["ref"]["species"]
build = config["ref"]["build"]
release = config["ref"]["release"]
genome_name = f"genome.{datatype_genome}.{species}.{build}.{release}"
genome_prefix = f"resources/{genome_name}"
genome = f"{genome_prefix}.fasta"
genome_fai = f"{genome}.fai"
genome_dict = f"{genome_prefix}.dict"
pangenome_name = f"pangenome.{species}.{build}"
pangenome_prefix = f"resources/{pangenome_name}"

# cram variables
use_cram = config.get("use_cram", False)
alignmend_ending = "cram" if use_cram else "bam"
alignmend_index_ending = "crai" if use_cram else "bai"
alignmend_ending_index_ending = "cram.crai" if use_cram else "bam.bai"

delly_excluded_regions = {
    ("homo_sapiens", "GRCh38"): "human.hg38",
    ("homo_sapiens", "GRCh37"): "human.hg19",
}


def _group_or_sample(row):
    group = row.get("group", None)
    if pd.isnull(group):
        return row["sample_name"]
    return group


samples["group"] = [_group_or_sample(row) for _, row in samples.iterrows()]

if "umi_read" not in samples.columns:
    samples["umi_read"] = pd.NA

validate(samples, schema="../schemas/samples.schema.yaml")


# Does this correctly return groups where fusion and variants are set?
def get_calling_groups(calling_type):
    return samples.loc[
        samples["calling"].str.contains(calling_type),
        "group",
    ].unique()


groups = samples["group"].unique()
calling_types = samples["calling"].str.split(",").explode().unique().tolist()
variants_groups = get_calling_groups("variants")
fusions_groups = get_calling_groups("fusions")

if "groups" in config:
    group_annotation = (
        pd.read_csv(config["groups"], sep="\t", dtype={"group": str})
        .set_index("group")
        .sort_index()
    )
    group_annotation = group_annotation.loc[groups]
else:
    group_annotation = pd.DataFrame({"group": groups}).set_index("group")

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


def is_activated(xpath):
    c = config
    for entry in xpath.split("/"):
        c = c.get(entry, {})
    return bool(c.get("activate", False))


custom_alignment_props = (
    (
        pd.read_csv(
            config["custom_alignment_properties"]["tsv"],
            sep="\t",
            dtype={"name": str, "path": str},
            comment="#",
        )
        .set_index(["name"], drop=False)
        .sort_index()
    )
    if is_activated("custom_alignment_properties")
    else None
)


def get_calling_events(calling_type):
    events = [
        event
        for event, entries in config["calling"]["fdr-control"]["events"].items()
        if calling_type in entries.get("types", ["variants"])
    ]
    return events


def get_heterogeneous_labels():
    nunique = group_annotation.nunique()
    cols_to_drop = nunique[nunique == 1].index
    return group_annotation.drop(cols_to_drop, axis=1).T


def get_final_output(wildcards):
    final_output = expand(
        "results/qc/multiqc/{group}.html",
        group=groups,
    )

    final_output.extend(
        expand(
            "results/datavzrd-report/{group}.coverage",
            group=groups,
        )
    )

    for calling_type in calling_types:
        if config["report"]["activate"]:
            for event in get_calling_events(calling_type):
                if lookup(
                    dpath=f"calling/fdr-control/events/{event}/report",
                    within=config,
                    default=True,
                ):
                    final_output.extend(
                        expand(
                            "results/datavzrd-report/{batch}.{event}.{calling_type}.fdr-controlled",
                            batch=get_report_batches(calling_type),
                            event=event,
                            calling_type=calling_type,
                        )
                    )
                else:
                    final_output.extend(
                        expand(
                            "results/final-calls/{group}.{event}.{calling_type}.fdr-controlled.bcf",
                            group=(
                                variants_groups
                                if calling_type == "variants"
                                else fusions_groups
                            ),
                            event=event,
                            calling_type=calling_type,
                        )
                    )
        else:
            final_output.extend(
                expand(
                    "results/final-calls/{group}.{event}.{calling_type}.fdr-controlled.bcf",
                    group=(
                        variants_groups
                        if calling_type == "variants"
                        else fusions_groups
                    ),
                    event=get_calling_events(calling_type),
                    calling_type=calling_type,
                )
            )

        if lookup(dpath="maf/activate", within=config, default=False):
            final_output.extend(
                expand(
                    "results/maf/{group}.{event}.{calling_type}.fdr-controlled.maf",
                    group=(
                        variants_groups
                        if calling_type == "variants"
                        else fusions_groups
                    ),
                    event=get_calling_events(calling_type),
                    calling_type=calling_type,
                )
            )

        if config["tables"]["activate"]:
            final_output.extend(
                expand(
                    "results/tables/{group}.{event}.{calling_type}.fdr-controlled.tsv",
                    group=(
                        variants_groups
                        if calling_type == "variants"
                        else fusions_groups
                    ),
                    event=get_calling_events(calling_type),
                    calling_type=calling_type,
                )
            )
            if config["tables"].get("generate_excel", False):
                final_output.extend(
                    expand(
                        "results/tables/{group}.{event}.{calling_type}.fdr-controlled.xlsx",
                        group=(
                            variants_groups
                            if calling_type == "variants"
                            else fusions_groups
                        ),
                        event=get_calling_events(calling_type),
                        calling_type=calling_type,
                    )
                )
    final_output.extend(get_mutational_burden_targets())
    final_output.extend(get_mutational_signature_targets())

    if is_activated("population/db"):
        final_output.append(lookup(dpath="population/db/path", within=config))

    return final_output


def get_gather_calls_input(ext="bcf"):
    def inner(wildcards):
        if wildcards.by == "odds":
            pattern = "results/calls/{{{{group}}}}.{{{{event}}}}.{{{{calling_type}}}}.{{scatteritem}}.filtered_odds.{ext}"
        elif wildcards.by == "ann":
            pattern = "results/calls/{{{{group}}}}.{{{{event}}}}.{{{{calling_type}}}}.{{scatteritem}}.filtered_ann.{ext}"
        else:
            raise ValueError(
                "Unexpected wildcard value for 'by': {}".format(wildcards.by)
            )
        return gather.calling(pattern.format(ext=ext))

    return inner


def get_control_fdr_input(wildcards):
    query = get_fdr_control_params(wildcards)
    if (
        not is_activated("benchmarking")
        and query["filter"]
        and wildcards.calling_type == "variants"
    ):
        by = "ann" if query["local"] else "odds"
        return "results/calls/{{group}}.{{event}}.{{calling_type}}.filtered_{by}.bcf".format(
            by=by
        )
    else:
        return "results/final-calls/{group}.{calling_type}.annotated.bcf"


def get_aligner(wildcards):
    if get_sample_datatype(wildcards.sample) == "rna":
        return "star"
    elif is_activated("ref/pangenome"):
        return "vg"
    else:
        return "bwa"


def get_fastp_input(wildcards):
    unit = units.loc[wildcards.sample].loc[wildcards.unit]

    if pd.isna(unit["fq1"]):
        return get_sra_reads(wildcards.sample, wildcards.unit, ["1", "2"])

    fq1 = get_raw_reads(unit.sample_name, unit.unit_name, "fq1")

    if len(fq1) == 1:

        def get_reads(fq):
            return get_raw_reads(unit.sample_name, unit.unit_name, fq)[0]

    else:
        ending = ".gz" if unit["fq1"].endswith("gz") else ""

        def get_reads(fq):
            return f"pipe/fastp/{unit.sample_name}/{unit.unit_name}.fq1.fastq{ending}"

    if pd.isna(unit["fq2"]):
        # single end sample
        return get_reads("fq1")
    else:
        # paired end sample
        return [get_reads("fq1"), get_reads("fq2")]


def get_sra_reads(sample, unit, fq):
    unit = units.loc[sample].loc[unit]
    # SRA sample (always paired-end for now)
    accession = unit["sra"]
    return expand("sra/{accession}_{read}.fastq.gz", accession=accession, read=fq)


def get_raw_reads(sample, unit, fq):
    pattern = units.loc[sample].loc[unit, fq]

    if pd.isna(pattern):
        assert fq.startswith("fq")
        fq = fq[len("fq") :]
        return get_sra_reads(sample, unit, fq)

    if type(pattern) is not str and len(pattern) > 1:
        raise ValueError(
            f"Multiple units.tsv entries found for sample '{sample}' and "
            f"unit '{unit}'.\n"
            "The units.tsv should contain only one entry for each combination "
            "of sample and unit.\n"
            "Found:\n"
            f"{pattern}"
        )

    if "*" in pattern:
        files = sorted(glob.glob(units.loc[sample].loc[unit, fq]))
        if not files:
            raise ValueError(
                "No raw fastq files found for unit pattern {} (sample {}). "
                "Please check your sample sheet.".format(unit, sample)
            )
    else:
        files = [pattern]
    return files


def get_fastp_pipe_input(wildcards):
    return get_raw_reads(wildcards.sample, wildcards.unit, wildcards.fq)


def get_fastqc_input(wildcards):
    return get_raw_reads(wildcards.sample, wildcards.unit, wildcards.fq)[0]


def get_fastp_adapters(wildcards):
    unit = units.loc[wildcards.sample].loc[wildcards.unit]
    try:
        adapters = unit["adapters"]
        if isinstance(adapters, str):
            # Autotrimming is enabled by default.
            # Therefore no adapter parameter needs to be passed.
            if adapters == "auto_trim":
                return ""
            else:
                return adapters
        return ""
    except KeyError:
        return ""


def get_fastp_extra(wildcards):
    extra = config["params"]["fastp"]
    if sample_has_umis(wildcards.sample):
        extra += get_annotate_umis_params(wildcards)
    return extra


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


def get_map_reads_input(wildcards):
    if is_paired_end(wildcards.sample):
        return [
            "results/merged/{sample}_R1.fastq.gz",
            "results/merged/{sample}_R2.fastq.gz",
        ]
    return "results/merged/{sample}_single.fastq.gz"


def get_star_reads_input(wildcards, r2=False):
    match (bool(is_paired_end(wildcards.sample)), r2):
        case (True, False):
            return "results/merged/{sample}_R1.fastq.gz"
        case (True, True):
            return "results/merged/{sample}_R2.fastq.gz"
        case (False, False):
            return "results/merged/{sample}_single.fastq.gz"
        case (False, True):
            return []


def get_group_aliases(group):
    return samples.loc[samples["group"] == group]["alias"]


def get_group_tumor_aliases(group):
    aliases = get_group_aliases(group)
    return aliases[aliases.str.startswith("tumor")]


def get_group_samples(group):
    return samples.loc[samples["group"] == group]["sample_name"]


def get_group_sample_aliases(wildcards, controls=True):
    if controls:
        return samples.loc[samples["group"] == wildcards.group]["alias"]
    return samples.loc[
        (samples["group"] == wildcards.group) & (samples["control"] == "no")
    ]["alias"]


def get_sample_datatype(sample):
    return samples.loc[[sample], "datatype"].iloc[0]


def get_markduplicates_input(wildcards):
    aligner = get_aligner(wildcards)
    if sample_has_umis(wildcards.sample):
        return f"results/mapped/{aligner}/{{sample}}.annotated.bam"
    else:
        return f"results/mapped/{aligner}/{{sample}}.sorted.bam"


def get_recalibrate_quality_input(wildcards, bai=False):
    ext = "bai" if bai else "bam"
    datatype = get_sample_datatype(wildcards.sample)
    if datatype == "rna":
        return "results/split/{{sample}}.{ext}".format(ext=ext)
    # Post-processing of DNA samples
    if is_activated("calc_consensus_reads"):
        return "results/consensus/{{sample}}.{ext}".format(ext=ext)
    else:
        return get_consensus_input(wildcards, bai)


def get_consensus_input(wildcards, bai=False):
    ext = "bai" if bai else "bam"
    if sample_has_primers(wildcards):
        return f"results/trimmed/{{sample}}.trimmed.{ext}"
    else:
        return get_trimming_input(wildcards, bai)


def get_trimming_input(wildcards, bai=False):
    ext = "bai" if bai else "bam"
    aligner = get_aligner(wildcards)
    if is_activated("remove_duplicates"):
        return "results/dedup/{{sample}}.{ext}".format(ext=ext)
    else:
        return "results/mapped/{aligner}/{{sample}}.sorted.{ext}".format(
            aligner=aligner, ext=ext
        )


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


def extract_unique_sample_column_value(sample, col_name):
    result = samples.loc[samples["sample_name"] == sample, col_name].drop_duplicates()
    if type(result) is not str:
        if len(result) > 1:
            ValueError(
                "If a sample is specified multiple times in a samples.tsv"
                "sheet, all columns except 'group' must contain identical"
                "entries across the occurrences (rows).\n"
                f"Here we have sample '{sample}' with multiple entries for"
                f"the '{col_name}' column, namely:\n"
                f"{result}\n"
            )
        else:
            result = result.squeeze()
    return result


def get_sample_primer_fastas(sample):
    if isinstance(primer_panels, pd.DataFrame):
        panel = extract_unique_sample_column_value(sample, "panel")
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
        panel = extract_unique_sample_column_value(wc.sample, "panel")
        return f"results/primers/{panel}_primer_regions.tsv"
    return "results/primers/uniform_primer_regions.tsv"


def get_markduplicates_extra(wc):
    c = config["params"]["picard"]["MarkDuplicates"]

    if sample_has_umis(wc.sample):
        b = "--BARCODE_TAG BX"
    else:
        b = ""

    if is_activated("calc_consensus_reads"):
        d = "--TAG_DUPLICATE_SET_MEMBERS true"
    else:
        d = ""

    return f"{c} {b} {d}"


def get_group_bams(wildcards, bai=False):
    ext = "bai" if bai else "bam"
    return expand(
        "results/recal/{sample}.{ext}",
        sample=get_group_samples(wildcards.group),
        ext=ext,
    )


def get_arriba_group_candidates(wildcards, csi=False):
    ext = ".csi" if csi else ""
    return expand(
        "results/candidate-calls/{sample}.arriba.bcf{ext}",
        sample=get_group_samples(wildcards.group),
        ext=ext,
    )


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


def get_all_group_observations(wildcards):
    return expand(
        "results/observations/{group}/{sample}.{caller}.all.bcf",
        caller=wildcards.caller,
        group=wildcards.group,
        sample=get_group_samples(wildcards.group),
    )


# TODO maybe use get_read_group with prefix="--outSAMattrRGline "?
def get_star_read_group(wildcards):
    """Denote sample name and platform in read group."""
    platform = extract_unique_sample_column_value(wildcards.sample, "platform")
    return r"--outSAMattrRGline ID:{sample} SM:{sample} PL:{platform}".format(
        sample=wildcards.sample, platform=platform
    )


def get_read_group(prefix: str):
    def inner(wildcards):
        """Denote sample name and platform in read group."""
        platform = extract_unique_sample_column_value(wildcards.sample, "platform")
        return r"{prefix}'@RG\tID:{sample}\tSM:{sample}\tPL:{platform}'".format(
            sample=wildcards.sample, platform=platform, prefix=prefix
        )

    return inner


def get_map_reads_sorting_params(wildcards, ordering=False):
    match (sample_has_umis(wildcards.sample), ordering):
        case (True, True):
            return "queryname"
        case (True, False):
            return "fgbio"
        case (False, True):
            return "coordinate"
        case (False, False):
            return "samtools"


def get_preprocessed_sorting_input(wildcards):
    if wildcards.aligner == "vg":
        return "results/mapped/vg/{sample}.rg.bam"
    else:
        return get_add_readgroup_input(wildcards)


def get_mutational_burden_targets():
    mutational_burden_targets = []
    if is_activated("mutational_burden"):
        for group in variants_groups:
            mutational_burden_targets.extend(
                expand(
                    "results/plots/mutational-burden/{group}.{alias}.{mode}.mutational-burden.svg",
                    group=group,
                    mode=config["mutational_burden"].get("mode", "curve"),
                    alias=get_group_tumor_aliases(group),
                )
            )
    return mutational_burden_targets


def get_mutational_signature_targets():
    mutational_signature_targets = []
    if is_activated("mutational_signatures"):
        for group in variants_groups:
            mutational_signature_targets.extend(
                expand(
                    "results/plots/mutational_signatures/{group}.{event}.html",
                    group=variants_groups,
                    event=config["mutational_signatures"].get("events"),
                )
            )
    return mutational_signature_targets


def get_scattered_calls(ext="bcf"):
    def inner(wildcards):
        caller = "arriba" if wildcards.calling_type == "fusions" else variant_caller
        return expand(
            "results/calls/{{group}}.{caller}.{{scatteritem}}.{ext}",
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


def get_annotated_bcf(wildcards, index=False):
    ext = ".csi" if index else ""
    selection = (
        get_selected_annotations() if wildcards.calling_type == "variants" else ""
    )
    return (
        "results/calls/{group}.{calling_type}.{scatteritem}{selection}.bcf{ext}".format(
            group=wildcards.group,
            calling_type=wildcards.calling_type,
            selection=selection,
            scatteritem=wildcards.scatteritem,
            ext=ext,
        )
    )


def get_gather_annotated_calls_input(ext="bcf"):
    def inner(wildcards):
        selection = (
            get_selected_annotations() if wildcards.calling_type == "variants" else ""
        )
        return gather.calling(
            "results/calls/{{{{group}}}}.{{{{calling_type}}}}.{{scatteritem}}{selection}.{ext}".format(
                ext=ext, selection=selection
            )
        )

    return inner


def get_candidate_calls(wc):
    filter = config["calling"]["filter"].get("candidates")
    if filter and wc.caller != "arriba":
        return "results/candidate-calls/{group}.{caller}.{scatteritem}.filtered.bcf"
    else:
        return "results/candidate-calls/{group}.{caller}.{scatteritem}.bcf"


def _get_report_batch(calling_type, batch):
    if batch == "all":
        _groups = variants_groups if calling_type == "variants" else fusions_groups
    else:
        _groups = samples.loc[
            samples[config["report"]["stratify"]["by-column"]] == batch,
            "group",
        ].unique()
    if not any(_groups):
        raise ValueError("No samples found. Is your sample sheet empty?")
    return _groups


def get_report_batch(calling_type):
    def inner(wildcards):
        return _get_report_batch(calling_type, wildcards.batch)

    return inner


def get_report_batches(calling_type):
    if is_activated("report/stratify"):
        yield "all"
        yield from samples[samples["calling"].str.contains(calling_type)][
            config["report"]["stratify"]["by-column"]
        ].unique()
    else:
        yield "all"


def get_merge_calls_input(ext="bcf"):
    def inner(wildcards):
        vartype = (
            ["SNV", "INS", "DEL", "MNV", "BND", "INV", "DUP", "REP"]
            if wildcards.calling_type == "variants"
            else ["BND"]
        )
        return expand(
            "results/calls/{{group}}.{vartype}.{{event}}.{{calling_type}}.fdr-controlled.{ext}",
            ext=ext,
            vartype=vartype,
        )

    return inner


def get_plugin_aux(plugin, index=False):
    if plugin in config["annotations"]["vep"]["final_calls"]["plugins"]:
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

    # local is a formerly supported key which we support as a fallback now
    # in order to a avoid a breaking change in the config.
    local = lookup(
        dpath="local",
        within=query,
        default=lookup(
            dpath="calling/fdr-control/local",
            within=config,
            default=True,
        ),
    )
    mode = lookup(
        dpath="mode",
        within=query,
        default=lookup(
            dpath="calling/fdr-control/mode",
            within=config,
            default="local-smart" if local else "global-smart",
        ),
    )

    mode = f"--mode {mode}"

    retain_artifacts = lookup(
        dpath="retain-artifacts",
        within=query,
        default=lookup(
            dpath="calling/fdr-control/retain-artifacts", within=config, default=False
        ),
    )

    retain_artifacts = "--smart-retain-artifacts" if retain_artifacts else ""

    return {
        "threshold": threshold,
        "events": events,
        "mode": mode,
        "retain_artifacts": retain_artifacts,
        "local": local,
        "filter": query.get("filter"),
    }


def get_fixed_candidate_calls(ext="bcf"):
    def inner(wildcards):
        if wildcards.caller == "delly":
            return expand(
                "results/candidate-calls/{{group}}.delly.no_bnds.{ext}",
                ext=ext,
            )
        else:
            return expand(
                "results/candidate-calls/{{group}}.{{caller}}.{ext}",
                ext=ext,
            )

    return inner


def get_filter_targets(wildcards, input):
    if input.predefined:
        return " | bedtools intersect -a /dev/stdin -b {input.predefined} ".format(
            input=input
        )
    else:
        return ""


def get_filter_expression(filter_name):
    filter_entry = config["calling"]["filter"][filter_name]
    if isinstance(filter_entry, str):
        return filter_entry
    else:
        return filter_entry["expression"]


def get_filter_aux_entries(filter_name):
    filter_entry = config["calling"]["filter"][filter_name]
    if isinstance(filter_entry, str):
        return {}
    else:
        aux = filter_entry.get("aux-files", {})
        return aux  # [f"--aux {name} {path}" for name, path in aux.items()]


def get_annotation_filter_names(wildcards):
    entry = config["calling"]["fdr-control"]["events"][wildcards.event]["filter"]
    filter_names = [entry] if isinstance(entry, str) else entry
    return filter_names


def get_annotation_filter_expression(wildcards):
    filters = [
        get_filter_expression(filter)
        for filter in get_annotation_filter_names(wildcards)
    ]
    return " and ".join(map("({})".format, filters)).replace('"', '\\"')


def get_annotation_filter_aux(wildcards, input):
    return [
        f"--aux {name}={path}"
        for filter in get_annotation_filter_names(wildcards)
        for name, path in zip(get_filter_aux_entries(filter).keys(), input.aux)
    ]


def get_annotation_filter_aux_files(wildcards):
    return [
        # aux files are considered as part of the config, hence local
        local(path)
        for filter_name in get_annotation_filter_names(wildcards)
        for name, path in get_filter_aux_entries(filter_name).items()
    ]


def get_candidate_filter_expression(wildcards):
    f = config["calling"]["filter"]["candidates"]
    if isinstance(f, dict):
        expression = f["expression"]
    else:
        expression = f
    return expression.replace('"', '\\"')


def get_candidate_filter_aux_files():
    if "candidates" not in config["calling"]["filter"]:
        return []
    else:
        # aux files are considered as part of the config, hence local
        return [local(path) for name, path in get_filter_aux_entries("candidates").items()]


def get_candidate_filter_aux(wildcards, input):
    if "candidates" not in config["calling"]["filter"]:
        return ""
    else:
        return [
            f"--aux {name}={path}"
            for name, path in zip(get_filter_aux_entries("candidates").keys(), input.aux)
        ]


def get_varlociraptor_obs_args(wildcards, input):
    return [
        "{}={}".format(s, f)
        for s, f in zip(get_group_aliases(wildcards.group), input.obs)
    ]


def get_varlociraptor_params(wildcards, params):
    if wildcards.caller == "arriba":
        params += " --propagate-info-fields GENE_NAME GENE_ID EXON_NUMBER"
    return params


wildcard_constraints:
    group="|".join(groups),
    sample="|".join(samples["sample_name"]),
    caller="|".join(["freebayes", "delly", "arriba"]),
    filter="|".join(config["calling"]["filter"]),
    event="|".join(config["calling"]["fdr-control"]["events"].keys()),
    regions_type="|".join(["expanded", "covered"]),
    calling_type="|".join(["fusions", "variants"]),


variant_caller = list(
    filter(
        None,
        [
            (
                "freebayes"
                if is_activated("calling/freebayes")
                and samples["calling"].str.contains("variants").any()
                else None
            ),
            (
                "delly"
                if is_activated("calling/delly")
                and samples["calling"].str.contains("variants").any()
                else None
            ),
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
                    f"SnpSift annotate -name '{prefix}_' {repr(path)} /dev/stdin"
                    for (prefix, _), path in zip(annotations, input.annotations)
                ]
            )
        )
    else:
        return ""


def get_cleaned_population_db(idx=False):
    if is_activated("population/db"):
        db = lookup(dpath="population/db/path", within=config)
        csi = ".csi" if idx else ""
        if not db.endswith(".bcf"):
            raise ValueError("Population database must be a BCF file.")
        if exists(db):
            return f"results/population_db.cleaned.bcf{csi}"
        else:
            return []
    return []


def get_population_bcfs(idx=False):
    csi = ".csi" if idx else ""
    return expand(
        "results/population/{group}.variants.filtered.bcf{idx}",
        group=variants_groups,
        idx=csi,
    )


def get_population_db(use_before_update=False):
    if is_activated("population/db"):
        db = lookup(dpath="population/db/path", within=config)
        if not db.endswith(".bcf"):
            raise ValueError("Population database must be a BCF file.")
        if use_before_update:
            if exists(db):
                return before_update(db)
            else:
                return []
        else:
            return update(db)
    return []


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


def get_untrimmed_fastqs(wc):
    return units.loc[units.sample_name == wc.sample, wc.read]


def get_trimmed_fastqs(wc):
    if units.loc[wc.sample, "adapters"].notna().all():
        return expand(
            "results/trimmed/{sample}/{unit}_{read}.fastq.gz",
            unit=units.loc[wc.sample, "unit_name"],
            sample=wc.sample,
            read=wc.read,
        )
    else:
        fq = "fq1" if wc.read == "R1" or wc.read == "single" else "fq2"
        return [
            read
            for unit in units.loc[wc.sample, "unit_name"]
            for read in get_raw_reads(wc.sample, unit, fq)
        ]


def get_annotation_fields_for_tables(wildcards):
    annotation_fields = [
        "Amino_acids",
        "CANONICAL",
        "CLIN_SIG",
        "Consequence",
        "Feature",
        "Gene",
        "gnomADg_AF",
        "HGVSc",
        "HGVSg",
        "HGVSp",
        "IMPACT",
        "MANE_PLUS_CLINICAL",
        "Protein_position",
        "SWISSPROT",
        "SYMBOL",
    ]

    annotation_fields.extend(
        [
            field
            for field in lookup(
                dpath="tables/output/annotation_fields", within=config, default=[]
            )
            if field not in annotation_fields
        ]
    )
    for plugin in ["REVEL", "SpliceAI", "AlphaMissense"]:
        if any(
            entry.startswith(plugin)
            for entry in config["annotations"]["vep"]["final_calls"]["plugins"]
        ):
            if plugin == "REVEL":
                annotation_fields.append("REVEL")
            elif plugin == "SpliceAI":
                annotation_fields.extend(
                    [
                        "SpliceAI_pred_DS_AG",
                        "SpliceAI_pred_DS_AL",
                        "SpliceAI_pred_DS_DG",
                        "SpliceAI_pred_DS_DL",
                    ]
                )
            elif plugin == "AlphaMissense":
                annotation_fields.append(
                    "am_pathogenicity",
                )

    return annotation_fields


def get_info_fusion_fields_for_tables(wildcards):
    if wildcards.calling_type == "fusions":
        return [
            "MATEID",
            "GENE_NAME",
            "GENE_ID",
            "EXON_NUMBER",
        ]
    else:
        return []


def get_format_fields_for_tables(wildcards):
    format_fields = ["AF", "DP"]

    if lookup(dpath="tables/output/short_observations", within=config, default=False):
        format_fields.extend(["SROBS", "SAOBS"])

    if lookup(dpath="tables/output/observations", within=config, default=False):
        format_fields.extend(
            [
                "OBS",
            ]
        )

    return format_fields


def get_info_prob_fields_for_tables(wildcards, input):
    if lookup(dpath="tables/output/event_prob", within=config, default=False):
        with open(input.scenario, "r") as scenario_file:
            scenario = yaml.load(scenario_file, Loader=yaml.SafeLoader)
            events = list(scenario["events"].keys())
            events += ["artifact", "absent"]
            return events
    else:
        return []


def get_vembrane_config(wildcards, input):
    # map VCF field names (keys) to wanted column names in the vembrane
    # output header (values); as we add them in batches by field type, and
    # not in the order we want them downstream, we explicitly order them
    # before this function returns the fields in the vembrane expression
    # and header format
    columns_dict = {
        "CHROM": "chromosome",
        "POS": "position",
        "REF": "reference allele",
        "ALT": "alternative allele",
        "INFO['END']": "end position",
        "INFO['EVENT']": "event",
        "ID": "id",
    }

    def append_items(items, renames, field_func, header_func=None):
        for item in items:
            if item in renames:
                header_name = renames[item]["name"]
                if "expr" in renames[item]:
                    parts_field = renames[item]["expr"]
                else:
                    parts_field = field_func(item)
            else:
                parts_field = field_func(item)
                header_name = header_func(item) if header_func else item
            columns_dict[parts_field] = header_name

    # FORMAT fields
    format_fields_names = {
        "AF": "allele frequency",
        "DP": "read depth",
        "SROBS": "short ref observations",
        "SAOBS": "short alt observations",
        "OBS": "observations",
    }

    format_fields = get_format_fields_for_tables(wildcards)
    aliases = get_group_sample_aliases(wildcards)

    for f in format_fields:
        append_items(
            aliases,
            dict(),
            f"FORMAT['{f}']['{{}}']".format,
            f"{{}}: {format_fields_names[f] if f in format_fields_names else f}".format,
        )

    # INFO fields
    rename_info_fields = {
        "MATEID": {"name": "mateid", "expr": "INFO['MATEID'][0]"},
        "GENE_NAME": {
            "name": "feature_name",
        },
        "GENE_ID": {
            "name": "feature_id",
        },
        "EXON_NUMBER": {"name": "exon", "expr": "INFO['EXON_NUMBER'][0]"},
    }

    ## INFO fields holding varlociraptor probabilities
    info_prob_fields = get_info_prob_fields_for_tables(wildcards, input)
    append_items(
        info_prob_fields,
        rename_info_fields,
        lambda x: f"INFO['PROB_{x.upper()}']",
        "prob: {}".format,
    )

    ## INFO fields relevant in fusion calling, only added for 'fusion' calling
    info_fusion_fields = get_info_fusion_fields_for_tables(wildcards)
    append_items(
        info_fusion_fields, rename_info_fields, "INFO['{}']".format, lambda x: x.lower()
    )

    ## INFO field annotations from the INFO['ANN'] field
    if wildcards.calling_type != "fusions":
        rename_ann_fields = {
            "Protein_position": {
                "name": "protein position",
                "expr": "ANN['Protein_position'].raw",
            },
            "Amino_acids": {
                "name": "protein alteration (short)",
            },
            "CLIN_SIG": {
                "name": "clinical significance",
            },
            "gnomADg_AF": {
                "name": "gnomad genome af",
            },
            "EXON": {
                "name": "exon",
                "expr": "ANN['EXON'].raw",
            },
            "SpliceAI_pred_DS_AG": {
                "name": "spliceai acceptor gain",
            },
            "SpliceAI_pred_DS_AL": {
                "name": "spliceai acceptor loss",
            },
            "SpliceAI_pred_DS_DG": {
                "name": "spliceai donor gain",
            },
            "SpliceAI_pred_DS_DL": {
                "name": "spliceai donor loss",
            },
            "am_pathogenicity": {
                "name": "alphamissense",
            },
        }
        annotation_fields = get_annotation_fields_for_tables(wildcards)
        append_items(
            annotation_fields,
            rename_ann_fields,
            "ANN['{}']".format,
            lambda x: x.lower(),
        )
    # determine a stable sort order, to avoid implicit assumptions about field
    # order; downstream scripts split-call-tables.py and
    # join_fusion_partner.py will remove columns they respectively don't need

    def get_and_sort_multi_entry_fields(field_prefix):
        fields = [k for k in columns_dict.keys() if k.startswith(field_prefix)]
        return sorted(fields)

    # MAIN variants columns
    sort_order = [
        # variants only
        "ANN['SYMBOL']",
        "ANN['IMPACT']",
        "ANN['HGVSp']",
        "ANN['HGVSc']",
        # variants only
        "ANN['Consequence']",
        "ANN['CLIN_SIG']",
        "ANN['gnomADg_AF']",
        "ANN['EXON'].raw",
        "ANN['REVEL']",
        # variants only, split-call-tables.py will select the column with the
        # highest score and will put it in the same place
        "ANN['SpliceAI_pred_DS_AG']",
        "ANN['SpliceAI_pred_DS_AL']",
        "ANN['SpliceAI_pred_DS_DG']",
        "ANN['SpliceAI_pred_DS_DL']",
        # variants only
        "ANN['am_pathogenicity']",
        # variants & fusions
        "CHROM",
        "POS",
        # variants only
        "REF",
        "ALT",
        "ANN['Protein_position'].raw",
        "ANN['Amino_acids']",
        "ANN['CANONICAL']",
        "ANN['MANE_PLUS_CLINICAL']",
        # fusions only
        "INFO['GENE_NAME']",
        "INFO['GENE_ID']",
    ]
    multi_entry_fields_prefix = [
        # fusions only
        "INFO['EXON_NUMBER'][0]",
        # MAIN column for fusions, COLLAPSED column for variants
        "INFO['PROB_",
        # MAIN column for fusions & variants
        # gets moved ahead of consequence column for variants
        "FORMAT['AF']",
        # MAIN column for fusions, COLLAPSED column for variants
        "FORMAT['DP']",
        # COLLAPSED columns
        # fusions only, join key in join_fusion_partner.py
        "INFO['MATEID'][0]",
        # variants & fusions
        "FORMAT['SROBS']",
        # variants & fusions
        "FORMAT['SAOBS']",
        # variants & fusions
        "FORMAT['OBS']",
    ]
    for p in multi_entry_fields_prefix:
        sort_order.extend(get_and_sort_multi_entry_fields(p))
    sort_order.extend(
        [
            # only needed for fusions, as a join key; dropped in
            # split-call-tables.py
            "ID",
            # only needed for indexing of variants in split-call-tables.py
            "ANN['Gene']",
            # only needed for variants in oncoprint.py
            "ANN['HGVSg']",
            # only needed for variants in datavzrd_variants_calls
            "ANN['Feature']",
            "INFO['END']",
        ]
    )
    # sort columns, keeping only those in sort_order
    sorted_columns_dict = {k: columns_dict[k] for k in sort_order if k in columns_dict}
    join_items = ", ".join
    return {
        "expr": join_items(sorted_columns_dict.keys()),
        "header": join_items(sorted_columns_dict.values()),
    }


def get_umi_fastq(wildcards):
    umi_read = extract_unique_sample_column_value(wildcards.sample, "umi_read")
    if umi_read in ["fq1", "fq2"]:
        return "results/untrimmed/{S}_{R}.sorted.fastq.gz".format(
            S=wildcards.sample, R=umi_read
        )
    elif umi_read == "both":
        return expand(
            "results/untrimmed/{S}_{R}.sorted.fastq.gz",
            S=wildcards.sample,
            R=["fq1", "fq2"],
        )


def sample_has_primers(wildcards):
    sample_name = wildcards.sample

    if config["primers"]["trimming"].get("primers_fa1") or (
        "panel" in samples.columns
        and samples.loc[samples["sample_name"] == sample_name, "panel"].notna().any()
    ):
        if not is_paired_end(sample_name):
            raise WorkflowError(
                f"Primer trimming is only available for paired-end data. Sample '{sample_name}' is not paired-end."
            )
        return True
    return False


def sample_has_umis(sample):
    return pd.notna(extract_unique_sample_column_value(sample, "umi_read"))


def get_annotate_umis_params(wildcards):
    translate_param = {"fq1": "read1", "fq2": "read2", "both": "per_read"}
    return "--umi --umi_loc {read} --umi_len {umi_len}".format(
        read=translate_param[
            extract_unique_sample_column_value(wildcards.sample, "umi_read")
        ],
        umi_len=str(extract_unique_sample_column_value(wildcards.sample, "umi_len")),
    )


def get_dgidb_datasources():
    if config["annotations"]["dgidb"].get("datasources", ""):
        return "-s {}".format(" ".join(config["annotations"]["dgidb"]["datasources"]))
    return ""


def get_filter_params(wc):
    if isinstance(get_panel_primer_input(wc.panel), list):
        return "-b -F 12"
    return "-b -F 4"


def get_single_primer_flag(wc):
    if not isinstance(get_sample_primer_fastas(wc.sample), list):
        return "--first-of-pair"
    return ""


def get_shortest_primer_length(primers):
    primers = primers if isinstance(primers, list) else [primers]
    # set to 32 to match bwa-mem default value considering offset of 2
    min_length = 32
    for primer_file in primers:
        with open(primer_file, "r") as p:
            min_primer = min(
                [len(p.strip()) for i, p in enumerate(p.readlines()) if i % 2 == 1]
            )
            min_length = min(min_length, min_primer)
    return min_length


def get_primer_extra(wc, input):
    extra = rf"-R '@RG\tID:{wc.panel}\tSM:{wc.panel}' -L 100"
    min_primer_len = get_shortest_primer_length(input.reads)
    # Check if shortest primer is below default values
    if min_primer_len < 32:
        extra += f" -T {min_primer_len - 2}"
    if min_primer_len < 19:
        extra += f" -k {min_primer_len}"
    return extra


def get_datavzrd_data(impact="coding"):
    calling_type = "variants"
    if impact == "fusions":
        impact = "fusions.joined"
        calling_type = "fusions"
    pattern = "results/tables/{group}.{event}.{impact}.fdr-controlled.tsv"

    def inner(wildcards):
        return expand(
            pattern,
            impact=impact,
            event=wildcards.event,
            group=get_report_batch(calling_type),
        )

    return inner


def get_oncoprint_input(wildcards):
    groups = get_report_batch("variants")
    return expand(
        "results/tables/{group}.{event}.coding.fdr-controlled.tsv",
        group=groups,
        event=wildcards.event,
    )


def get_variant_oncoprint_tables(wildcards, input):
    if input.variant_oncoprints:
        oncoprint_dir = input.variant_oncoprints
        valid = re.compile(r"^[^/]+\.tsv$")
        tables = [f for f in os.listdir(oncoprint_dir) if valid.match(f)]
        assert all(table.endswith(".tsv") for table in tables)
        genes = [gene_table[:-4] for gene_table in tables]
        return list(
            zip(genes, expand(f"{oncoprint_dir}/{{oncoprint}}", oncoprint=tables))
        )
    else:
        return []


def get_datavzrd_report_labels(wildcards):
    event = config["calling"]["fdr-control"]["events"][wildcards.event]
    labels = {"batch": wildcards.batch}
    if "labels" in event:
        labels.update({key: str(value) for key, value in event["labels"].items()})
    else:
        labels["callset"] = wildcards.event.replace("_", " ")
    return labels


def get_datavzrd_report_subcategory(wildcards):
    event = config["calling"]["fdr-control"]["events"][wildcards.event]
    return event.get("subcategory", None)


def get_fastqc_results(wildcards):
    group_samples = get_group_samples(wildcards.group)
    sample_units = units.loc[group_samples]
    sra_units = pd.isna(sample_units["fq1"])
    paired_end_units = sra_units | ~pd.isna(sample_units["fq2"])

    # fastqc
    pattern = "results/qc/fastqc/{unit.sample_name}/{unit.unit_name}.{fq}_fastqc.zip"
    yield from expand(pattern, unit=sample_units.itertuples(), fq="fq1")
    yield from expand(
        pattern, unit=sample_units[paired_end_units].itertuples(), fq="fq2"
    )

    # fastp
    if sample_units["adapters"].notna().all():
        pattern = "results/trimmed/{unit.sample_name}/{unit.unit_name}.{mode}.qc.html"
        yield from expand(
            pattern, unit=sample_units[paired_end_units].itertuples(), mode="paired"
        )
        yield from expand(
            pattern, unit=sample_units[~paired_end_units].itertuples(), mode="single"
        )

    # samtools idxstats
    yield from expand(
        "results/qc/{sample}.bam.idxstats",
        sample=group_samples,
    )

    # samtools stats
    yield from expand(
        "results/qc/{sample}.bam.stats",
        sample=group_samples,
    )


def get_variant_oncoprints(wildcards):
    if len(_get_report_batch("variants", wildcards.batch)) > 1:
        return "results/tables/oncoprints/{wildcards.batch}.{wildcards.event}/variant-oncoprints"
    else:
        return []


def get_oncoprint(oncoprint_type):
    def inner(wildcards):
        if len(_get_report_batch("variants", wildcards.batch)) > 1:
            oncoprint_path = (
                f"results/tables/oncoprints/{wildcards.batch}.{wildcards.event}"
            )
            if oncoprint_type == "gene":
                return f"{oncoprint_path}/gene-oncoprint.tsv"
            elif oncoprint_type == "variant":
                return f"{oncoprint_path}/variant-oncoprints"
            else:
                raise ValueError(f"bug: unsupported oncoprint type {oncoprint_type}")
        else:
            return []

    return inner


def get_delly_excluded_regions():
    custom_excluded_regions = config["calling"]["delly"].get("exclude_regions", "")
    if custom_excluded_regions:
        return custom_excluded_regions
    elif delly_excluded_regions.get((species, build), False):
        return "results/regions/{species_build}.delly_excluded.bed".format(
            species_build=delly_excluded_regions[(species, build)]
        )
    else:
        return []


def get_alignment_props(wildcards):
    if is_activated("custom_alignment_properties"):
        alignment_prop_column = config["custom_alignment_properties"]["column"]
        prop_name = extract_unique_sample_column_value(
            wildcards.sample, alignment_prop_column
        )
        if pd.notna(prop_name):
            return custom_alignment_props.loc[prop_name, "path"]
    return f"results/alignment-properties/{wildcards.group}/{wildcards.sample}.json"


def get_pangenome_url(datatype):
    build = config["ref"]["build"].lower()
    source = config["ref"]["pangenome"]["source"]
    version = config["ref"]["pangenome"]["version"]
    if config["ref"]["species"] != "homo_sapiens" or build not in ["grch37", "grch38"]:
        raise ValueError(
            "Unsupported combination of species and build. Only homo_sapiens and GRCh37/GRCh38 are supported for pangenome mapping."
        )
    if source != "hprc":
        raise ValueError(
            "Unsupported pangenome source. Only 'hprc' is currently supported."
        )
    return f"https://s3-us-west-2.amazonaws.com/human-pangenomics/pangenomes/freeze/freeze1/minigraph-cactus/hprc-{version}-mc-{build}/hprc-{version}-mc-{build}.{datatype}"
