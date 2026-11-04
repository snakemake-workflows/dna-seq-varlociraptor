"""Check the candidate filter of config/config.yaml on annotated calls.

The .test configs define their own filters, hence the workflow tests do not
exercise the defaults that users start from.
"""

from pathlib import Path
import subprocess
import sys

import yaml

HERE = Path(__file__).parent
CALLS = HERE / "calls.vcf"
GENES = HERE / "genes.tsv"

MUST_BE_KEPT = {"high_modifier", "moderate", "high_moderate"}
EXPECTED_REMOVED = {"low", "modifier", "low_modifier"}


def get_expression(entry):
    """Return the vembrane expression of a filter entry of the config."""
    return entry if isinstance(entry, str) else entry["expression"]


def get_aux_names(entry):
    """Return the names of the aux files of a filter entry of the config."""
    return [] if isinstance(entry, str) else entry.get("aux-files", {}).keys()


def run_filter(expression, aux_names=()):
    """Apply the expression to the fixture and return the IDs of the kept records."""
    # The gene list of the config is an example, use one that matches the fixture.
    aux = [arg for name in aux_names for arg in ("--aux", f"{name}={GENES}")]
    filtered = subprocess.run(
        ["vembrane", "filter", *aux, expression, CALLS],
        check=True,
        stdout=subprocess.PIPE,
        text=True,
    )
    return {
        line.split("\t")[2]
        for line in filtered.stdout.splitlines()
        if not line.startswith("#")
    }


config = yaml.safe_load((HERE.parent.parent / "config/config.yaml").read_text())
filters = config["calling"]["filter"]
candidates = run_filter(
    get_expression(filters["candidates"]), get_aux_names(filters["candidates"])
)

errors = []
if MUST_BE_KEPT - candidates:
    errors.append(
        f"candidates removes HIGH or MODERATE variants: {sorted(MUST_BE_KEPT - candidates)}"
    )
if EXPECTED_REMOVED & candidates:
    errors.append(
        "candidates keeps variants with only LOW or MODIFIER annotations: "
        f"{sorted(EXPECTED_REMOVED & candidates)}. "
        "If this is intended, update EXPECTED_REMOVED."
    )

for name, event in config["calling"]["fdr-control"]["events"].items():
    entry = event.get("filter", [])
    names = [entry] if isinstance(entry, str) else entry
    if not names:
        continue
    kept = run_filter(
        " and ".join(f"({get_expression(filters[f])})" for f in names),
        [aux for f in names for aux in get_aux_names(filters[f])],
    )
    if not kept:
        errors.append(f"event {name} keeps no variant of the fixture")
    if kept - candidates:
        errors.append(
            f"candidates removes variants that event {name} accepts: {sorted(kept - candidates)}"
        )

if errors:
    sys.exit("\n".join(errors))
print("OK:", sorted(candidates))
