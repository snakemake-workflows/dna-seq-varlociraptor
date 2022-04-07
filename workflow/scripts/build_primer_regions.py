import pandas as pd


def parse_bed(log_file, out):
    print("chrom\tleft_start\tleft_end\tright_start\tright_end", file=out)
    for data_primers in pd.read_csv(
        snakemake.input[0],
        sep="\t",
        header=None,
        chunksize=chunksize,
        usecols=[0, 1, 2, 5],
    ):
        for row in data_primers.iterrows():
            row_id = row[0]
            row = row[1]
            if row[5] == "+":
                print(
                    "{chrom}\t{start}\t{end}\t-1\t-1".format(
                        chrom=row[0], start=row[1] + 1, end=row[2]
                    ),
                    file=out,
                )
            elif row[5] == "-":
                print(
                    "{chrom}\t-1\t-1\t{start}\t{end}".format(
                        chrom=row[0], start=row[1] + 1, end=row[2]
                    ),
                    file=out,
                )
            else:
                print("Invalid strand in row {}".format(row_id), file=log_file)


def parse_bedpe(log_file, out):
    for data_primers in pd.read_csv(
        snakemake.input[0],
        sep="\t",
        header=None,
        chunksize=chunksize,
        usecols=[0, 1, 2, 3, 4, 5],
    ):
        valid_primers = data_primers[0] == data_primers[3]
        valid_data = data_primers[valid_primers].copy()
        valid_data.iloc[:, [1, 4]] += 1
        valid_data.drop(columns=[3], inplace=True)
        valid_data.dropna(how="all", inplace=True)
        valid_data.to_csv(
            out,
            sep="\t",
            index=False,
            header=["chrom", "left_start", "left_end", "right_start", "right_end"],
        )
        print(
            data_primers[~valid_primers].to_csv(sep="\t", index=False, header=False),
            file=log_file,
        )


chunksize = 10**6
with open(snakemake.output[0], "w") as out:
    with open(snakemake.log[0], "w") as log_file:
        if snakemake.input[0].endswith("bedpe"):
            parse_bedpe(log_file, out)
        else:
            parse_bed(log_file, out)
