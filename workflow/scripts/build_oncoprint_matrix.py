import pysam
import pandas as pd
import os


input_files = snakemake.input

df = pd.DataFrame(columns=["Sample"])
for sample_file in input_files:
    variant_file = pysam.VariantFile(sample_file)
    sample_name = os.path.basename(sample_file).split(".")[0]
    gene_variant_dict = {"Sample": [sample_name]}
    for rec in variant_file.fetch():
        for sample in rec.samples:
            allele_frequencies = rec.samples[sample]["AF"] #Can be multiple entries
            for allele_frequency in allele_frequencies:
                variant = rec.info["SVLEN"]
                if variant[0]:
                    variant_type = "INDEL"
                else:
                    variant_type = "SNV" 
                transcripts = rec.info["ANN"]
                for transcript in transcripts:
                    gene = transcript.split("|")[3]
                    if gene not in gene_variant_dict:
                        gene_variant_dict[gene] = set()
                    gene_variant_dict[gene].add(variant_type)
                break
    for key, value in gene_variant_dict.items():
        gene_variant_dict[key] = ','.join(value)
    sample_df = pd.DataFrame(gene_variant_dict, index=[0])
    df = pd.concat([df, sample_df], join="outer", ignore_index=False, sort=False)
df.set_index("Sample", inplace=True)
with open(snakemake.output[0], 'w') as output_f:
    print(df.to_csv(sep="\t", index=True), file=output_f)
