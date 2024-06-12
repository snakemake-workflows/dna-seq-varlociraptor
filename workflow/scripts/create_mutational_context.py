import sys
import pysam
from Bio import SeqIO


def get_ref_triplet(ref_seq, variant_pos):
    left_idx = variant_pos - 1 if variant_pos > 0 else 0
    right_idx = variant_pos + 2
    return ref_seq[left_idx:right_idx]


sys.stderr = open(snakemake.log[0], "w")

output = open(snakemake.output[0], "w")
max_vaf = int(snakemake.wildcards.vaf) / 100

for bcf_path in snakemake.input.bcfs:
    reference = SeqIO.parse(snakemake.input.ref, "fasta")
    sample_name = bcf_path.split(".")[0]
    bcf = pysam.VariantFile(bcf_path)
    current_chrom_id = None
    current_chrom_seq = None
    for bcf_record in bcf:
        if float(bcf_record.samples["tumor"]["AF"][0]) > max_vaf:
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
        print(f"{ref_triplet}\t{alt_bases[0]}\t{sample_name}", file=output)
    output.close()
