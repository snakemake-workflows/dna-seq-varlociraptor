import pandas as pd


def get_vartype(rec):
    ref_allele = rec["reference allele"]
    alt_allele = rec["alternative allele"]
    if alt_allele == "<DEL>" or (len(ref_allele) > 1 and len(alt_allele) == 1):
        return "deletion"
    elif alt_allele == "<INS>" or (len(alt_allele) > 1 and len(ref_allele) == 1):
        return "insertion"
    elif alt_allele == "<INV>":
        return "inversion"
    elif alt_allele == "<DUP>":
        return "duplication"
    elif alt_allele == "<TDUP>":
        return "tandem duplication"
    elif alt_allele == "<CNV>":
        return "copy number variation"
    elif alt_allele.startswith("<"):
        # catch other special alleles
        return "complex"
    elif len(alt_allele) == 1 and len(ref_allele) == 1:
        return "snv"
    elif len(alt_allele) == len(ref_allele):
        return "mnv"
    else:
        return "breakend"


def load_calls(path, group):
    calls = pd.read_csv(
        path, sep="\t", usecols=["gene", "reference allele", "alternative allele"]
    )
    calls["vartype"] = calls.apply(get_vartype, axis="columns")
    calls["group"] = group
    return calls[["group", "gene", "vartype"]].drop_duplicates()


pd.concat(
    [
        load_calls(path, sample)
        for path, sample in zip(snakemake.input, snakemake.params.groups)
    ]
).to_csv(snakemake.output[0], sep="\t", index=False)
