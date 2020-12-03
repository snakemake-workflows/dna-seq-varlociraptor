#!/usr/bin/env python

import sys
import pandas as pd

infile = sys.argv[1]
outfile = sys.argv[2]

d = pd.read_csv(infile, sep="\t")

## print(d)

## coverage
## d["FATHERS_N_COV"] = d["FATHERS_N_REF"] + d["FATHERS_N_ALT"]
## d["MOTHERS_N_COV"] = d["MOTHERS_N_REF"] + d["MOTHERS_N_ALT"]

## variant allele frequency (VAF)
d["FATHERS_VAF"] = d["FATHERS_N_VAR"] / d["FATHERS_N_COV"]
d["MOTHERS_VAF"] = d["MOTHERS_N_VAR"] / d["MOTHERS_N_COV"]

d = d.fillna("NaN")

## reorder columns
new_column_order = ["ID",
    "FATHERS_N_COV", "FATHERS_N_VAR", "FATHERS_VAF", 
    "MOTHERS_N_COV",  "MOTHERS_N_VAR", "MOTHERS_VAF", 
    "PROB_FATHERS_ONLY", "PROB_MOTHERS_ONLY", "PROB_FATHERS_AND_MOTHERS"]
d = d.reindex(new_column_order, axis=1)

d.to_csv(outfile, sep="\t", index=False)

