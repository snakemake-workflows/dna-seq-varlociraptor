import os
import shutil
import sys
import tempfile

sys.stderr = open(snakemake.log[0], "w")

import subprocess as sp
import pysam

db_path = snakemake.output[0]


def sample_in_db(group):
    samples = pysam.VariantFile(db_path).header.samples
    if group in samples:
        if len(samples) == 1:
            return "only"
        else:
            return "one"
    else:
        return "not"


def copy_to_db(bcf_path):
    shutil.copyfile(bcf_path, db_path)


def update_db(bcf_path, group):
    sp.run(["bcftools", "index", "-f", db_path])

    match sample_in_db(group):
        case "only":
            copy_to_db(bcf_path)
            return
        case "one":
            # This should probably remove the existing sample
            _, db_tmp = tempfile.mkstemp()
            sp.run(
                [
                    "bcftools",
                    "view",
                    db_path,
                    "--samples",
                    f"^{group}",
                    "-Ob",
                    "-o",
                    db_tmp,
                ]
            )
            sp.run(["bcftools", "index", db_tmp])
        case "not":
            db_tmp = db_path

    sp.run(
        [
            "bcftools",
            "merge",
            "-m",
            "none",
            db_tmp,
            bcf_path,
            "-Ob",
            "-o",
            f"{db_path}.updated",
        ],
        check=True,
        stderr=sp.PIPE,
    )
    shutil.move(f"{db_path}.updated", db_path)


for bcf_path in snakemake.input.bcf:
    if os.path.exists(db_path):
        group = bcf_path.split("/")[2].split(".")[0]
        update_db(bcf_path, group)
    else:
        copy_to_db(bcf_path)
