import os
import argparse
from multiprocessing import Pool
from pysam import VariantFile
import tempfile
import sys

parser = argparse.ArgumentParser()
parser.add_argument("file_vcf")
parser.add_argument("command")
parser.add_argument("--stepsize", "-s", type=int, default=1000000)
parser.add_argument("--threads", "-t", type=int, default=8)
parser.add_argument("--tmp", default=None, help="Define the base directory for the tmp files")
args = parser.parse_args()
tmp_base = args.tmp
stepsize = args.stepsize

# create regions
f = VariantFile(args.file_vcf)
regions = []
for name, contig in f.header.contigs.items():
    last = 0
    for x in list(range(0, contig.length, stepsize)):
        stop = min(last + stepsize - 1, contig.length)
        regions.append(f"{name}:{last}-{stop}")
        last = stop+1
f.close() 

#command = "SnpSift annotate {ss_args} {args.database_vcf} /dev/stdin"
command = args.command

with tempfile.TemporaryDirectory(dir=tmp_base) as t:
    def f(x):
        #print("start", x, file=sys.stderr)
        cmd = f"bcftools view -r {x} {args.file_vcf} | bcftools view -t {x} | {command} | bcftools view -Ob > {t}/{x}.bcf"
        os.system(cmd)
        #print("done", x, file=sys.stderr)
        return f"{t}/{x}.bcf"

    with Pool(args.threads) as p:
        first = True
        for filename in p.imap(f, regions):
            print("merge", filename, file=sys.stderr)
            if first:
                cmd = f"bcftools view {filename}"
            else:
                cmd = f'bcftools view {filename} | grep -v "#"'
            os.system(cmd)
            os.remove(filename)
            first = False
