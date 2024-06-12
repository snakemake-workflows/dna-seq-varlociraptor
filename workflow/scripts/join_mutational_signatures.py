import pandas as pd

sys.stderr = open(snakemake.log[0], "w")

dataframes = []

for file_path in snakemake.input:
    max_vaf = float(filepath.rsplit(".", 2)[-2])
    df = pd.read_csv(file_path, sep="\t")
    df["max_vaf"] = max_vaf / 100
    dataframes.append(df)
final_signatures_df = pd.concat(dataframes, ignore_index=True)
final_signatures_df.to_csv(snakemake.output, sep="\t", index=False)
