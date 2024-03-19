import os
import shutil
import sys
import tempfile

sys.stderr = open(snakemake.log[0], "w")

import subprocess as sp
import pysam

db_path = snakemake.output[0]
db_path_cleaned = snakemake.input.cleaned_db


def db_empty():
    is_empty = True
    if os.path.exists(db_path):
        sp.run(["bcftools", "index", "-f", db_path])
        samples = pysam.VariantFile(db_path).header.samples
        if len(samples) != 0:
            is_empty = False
    return is_empty


def copy_to_db(bcf_path):
    shutil.copyfile(bcf_path, db_path)


def update_db(bcf_path):
    sp.run(["bcftools", "index", "-f", db_path])

    sp.run(
        [
            "bcftools",
            "merge",
            "-m",
            "none",
            db_path,
            bcf_path,
            "-Ob",
            "-o",
            f"{db_path}.updated",
        ],
        check=True,
        stderr=sp.PIPE,
    )
    shutil.move(f"{db_path}.updated", db_path)


if db_path_cleaned:
    copy_to_db(db_path_cleaned)

for bcf_path in snakemake.input.bcf:
    if db_empty():
        copy_to_db(bcf_path)
    else:
        update_db(bcf_path)
