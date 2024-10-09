#!env bash

# Replace convert_fusions_to_vcf.sh script by updated one until new version of arriba is released
curl https://raw.githubusercontent.com/suhrig/arriba/6fff196/scripts/convert_fusions_to_vcf.sh > $CONDA_PREFIX/bin/convert_fusions_to_vcf.sh
