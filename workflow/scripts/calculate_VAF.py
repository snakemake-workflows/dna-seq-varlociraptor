#!/usr/bin/env python

import sys
import pandas as pd

infile = sys.argv[1]
outfile = sys.argv[2]

d = pd.read_csv(infile, sep="\t")

print(d)

## coverage
# d["FATHER_N_COV"] = d["FATHER_N_REF"] + d["FATHER_N_ALT"]
# d["MOTHER_N_COV"] = d["MOTHER_N_REF"] + d["MOTHER_N_ALT"]

## variant allele frequency (VAF)
d["FATHER_VAF"] = d["FATHER_N_VAR"] / d["FATHER_N_COV"]
d["MOTHER_VAF"] = d["MOTHER_N_VAR"] / d["MOTHER_N_COV"]

d = d.fillna("NA")

## reorder columns
new_column_order = ["ID",
    "FATHER_N_COV", "FATHER_N_VAR", "FATHER_VAF", 
    "MOTHER_N_COV",  "MOTHER_N_VAR", "MOTHER_VAF", 
    "PROB_FATHER_ONLY", "PROB_MOTHER_ONLY", "PROB_FATHER_AND_MOTHER"]
d = d.reindex(new_column_order, axis=1)

d.to_csv(outfile, sep="\t", index=False)

