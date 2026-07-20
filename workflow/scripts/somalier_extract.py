import os
import tempfile
from snakemake.shell import shell

with tempfile.NamedTemporaryFile(dir="/dev/shm", suffix=".bam", delete=False) as microbam:
    microbam_path = microbam.name

try:
    # use samtools to constrain input bam, it can read faster on non SSD (network, HDD)
    shell(
        "samtools view -u -M -L "
        "<(zcat {snakemake.input.sites} | awk '!/^#/{{print $1\"\\t\"$2-1\"\\t\"$2}}') "
        "-o {microbam_path} {snakemake.input.bam} 2> {snakemake.log}"
    )
    shell("samtools index {microbam_path} 2>> {log}")
    shell(
        "somalier extract -d {snakemake.params.outdir} "
        "--sites {snakemake.input.sites} --sample-prefix {snakemake.wildcards.sample} "
        "-f {snakemake.input.ref} {microbam_path} 2>> {log}"
    )
finally:
    os.remove(microbam_path)