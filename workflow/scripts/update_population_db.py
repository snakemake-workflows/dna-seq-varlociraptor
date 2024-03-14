import os
import shutil
import sys
sys.stderr = open(snakemake.log[0], "w")

from pathlib import Path
import subprocess as sp
import pysam

db_path = snakemake.output[0]


def sample_in_db():
    samples = pysam.VariantFile(db_path).header.samples
    if snakemake.wildcards.group in samples:
        if len(samples) == 1:
            return "only"
        else:
            return "one"
    else:
        return "not"


def copy_to_db():
    shutil.copyfile(snakemake.input[0], db_path)


def update_db():
    match sample_in_db():
        case "only":
            copy_to_db()
            return
        case "one":
            aux_args = ["--samples", f"^{snakemake.wildcards.group}"]
        case "not":
            aux_args = []
    sp.run(["bcftools", "merge", "-m", "none", db_path, snakemake.input[0], f"{db_path}.updated"] + aux_args, check=True, stderr=sp.PIPE)
    shutil.move(f"{db_path}.updated", db_path)


if os.path.exists(db_path):
    update_db()
else:
    copy_to_db()