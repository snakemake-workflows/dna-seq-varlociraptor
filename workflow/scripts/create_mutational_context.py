import sys
import pysam
import pandas as pd
from Bio import SeqIO

sys.stderr = open(snakemake.log[0], "w")

def get_ref_triplet(ref_seq, variant_pos):
    left_idx = variant_pos - 1 if variant_pos > 0 else 0
    right_idx = variant_pos + 2
    return ref_seq[left_idx:right_idx]

min_vaf = int(snakemake.wildcards.vaf) / 100

bcf_path = snakemake.input.bcf
reference = SeqIO.parse(snakemake.input.ref, "fasta")
sample_name = bcf_path.split("/")[-1].split(".")[0]
bcf = pysam.VariantFile(bcf_path)
current_chrom_id = None
current_chrom_seq = None
single_base_substitutions = []

for bcf_record in bcf:
    if float(bcf_record.samples["tumor"]["AF"][0]) < min_vaf:
        continue
    variant_chrom = bcf_record.chrom
    variant_pos = bcf_record.pos - 1
    ref_base = bcf_record.ref
    alt_bases = bcf_record.alts
    # alternatively check ANN field for VARIANT_CLASS == SNV?
    if len(alt_bases) > 1:
        print("Record has mutliple alterations", file=sys.stderr)
        print(f"{variant_chrom}\t{variant_pos}", file=sys.stderr)
        exit()
    if len(ref_base) != 1 or len(alt_bases[0]) != 1:
        print(
            f"Record skipped - No SNV: {variant_chrom}\t{variant_pos}",
            file=sys.stderr,
        )
        continue
    while variant_chrom != current_chrom_id:
        current_chrom = next(reference)
        current_chrom_id = current_chrom.id
        current_chrom_seq = current_chrom.seq
    ref_triplet = get_ref_triplet(current_chrom_seq, variant_pos)
    if ref_triplet[1] != ref_base:
        print(
            f"Error: Missmatching reference base: {variant_chrom}\t{variant_pos}\t{ref_base}\t{ref_triplet}",
            file=sys.stderr,
        )
        exit()
    single_base_substitutions.append((ref_triplet, alt_bases[0]))

df = pd.DataFrame(single_base_substitutions)
df[2] = sample_name
mutation_count = len(df.index)
df.to_csv(snakemake.output.context, sep="\t", index=False, header=False)
with open(snakemake.output.counts, "w") as count_output:
    print(f"{min_vaf}\t{mutation_count}", file=count_output)