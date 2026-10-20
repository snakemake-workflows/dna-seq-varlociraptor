import json
import os
import sys
from urllib.request import urlopen

sys.stderr = open(snakemake.log[0], "w")

species = snakemake.params.species
url = f"https://rest.ensembl.org/info/genomes/{species}?content-type=application/json"
with urlopen(url) as response:
    taxid = json.load(response)["taxonomy_id"]

# input file format of ganon build-custom: file, target, taxid
with open(snakemake.output[0], "w") as out:
    print(os.path.abspath(snakemake.input.genome), species, taxid, sep="\t", file=out)
