import sys
import pysam
import pandas as pd

sys.stderr = open(snakemake.log[0], "w")

bcf_file = pysam.VariantFile(snakemake.input[0])
group = snakemake.wildcards.group

# Lists to store variant data
variants_data = []

for record in bcf_file:
    ref = record.ref
    alts = record.alts
    if len(ref) == 1 and all(len(alt) == 1 for alt in alts):
        chrom = record.chrom
        pos = record.pos
        # Construct a tuple representing the variant
        variant = (group, chrom, ref, alts[0])
        # Append variant data to the list
        variants_data.append(variant)

# Create a DataFrame from the collected variant data
variants_df = pd.DataFrame(variants_data, columns=["group", "chrom", "ref", "alt"])

# Drop duplicates based on 'chrom', 'ref', and 'alt'
variants_df.drop_duplicates(subset=["chrom", "ref", "alt"], inplace=True)

# Write the DataFrame to a TSV file
variants_df.to_csv(snakemake.output[0], sep="\t", index=False)