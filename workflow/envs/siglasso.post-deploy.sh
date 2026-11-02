#!env bash
set -o pipefail
Rscript --no-environ -e 'Sys.setenv(R_REMOTES_NO_ERRORS_FROM_WARNINGS="false"); remotes::install_github("gersteinlab/siglasso", ref = "9c5e6b07dcf0f9ebf1303c4c541c3dac9267de71", upgrade = "never")'
