import pandas as pd

chunksize = 10 ** 6
with open(snakemake.output[0], "w") as out:
    with open(snakemake.log[0], "w") as log_file:
        for data_primers in pd.read_csv(
            snakemake.input[0],
            sep="\t",
            header=None,
            chunksize=chunksize,
            usecols=[0, 1, 3, 5, 6],
        ):
            valid_primers = data_primers[0] == data_primers[3]
            print(
                data_primers[valid_primers]
                .drop(columns=[3, 6])
                .to_csv(sep="\t", index=False, header=False),
                file=out,
            )
            print(
                data_primers[~valid_primers].to_csv(
                    sep="\t", index=False, header=False
                ),
                file=log_file,
            )
