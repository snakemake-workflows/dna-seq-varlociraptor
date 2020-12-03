#!/usr/bin/env python

import sys
import pandas as pd

infile_upload = sys.argv[1]
infile_vaf = sys.argv[2]
outfile = sys.argv[3]

du = pd.read_csv(infile_upload, sep="\t")

## add ID column for merge
du["ID"] = du[["Chr", "Start", "End", "Ref", "Alt"]].astype(str).apply(lambda x: '_'.join(x), axis=1)
## print(du)

dv = pd.read_csv(infile_vaf, sep="\t")
## print(dv)

dm = pd.merge(du, dv, on="ID")
dm = dm.drop(columns=["ID"]) ## remove ID column
## print(dm)

## print(dm.to_string(index=False))
dm.to_csv(outfile, sep="\t", index=False, na_rep="NA")

