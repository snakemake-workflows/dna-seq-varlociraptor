import pandas as pd

sys.stderr = open(snakemake.log[0], "w")

dataframes = []

for file_path in snakemake.input:
    min_vaf = float(file_path.rsplit(".", 2)[-2])
    df = pd.read_csv(file_path, sep="\t")
    df.rename(columns={snakemake.wildcards.group: "Frequency"}, inplace=True)
    df["Frequency"] = df["Frequency"].astype(float)
    df["min_vaf"] = min_vaf / 100
    df = df[df["Frequency"] != 0]
    dataframes.append(df)
final_signatures_df = pd.concat(dataframes, ignore_index=True)
final_signatures_df.to_csv(snakemake.output[0], sep="\t", index=False)
