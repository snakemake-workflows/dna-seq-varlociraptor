import SigProfilerAssignment as spa
from SigProfilerAssignment import Analyzer as Analyze
# from SigProfilerMatrixGenerator import install as genInstall

# Install reference required?
# genInstall.install(snakemake.params.build)

Analyze.cosmic_fit(samples=snakemake.input[0], #maybe only folder supported?!
                   output=snakemake.output[0],
                   input_type="vcf", # Probably no bcf supported
                   genome_build=snakemake.params.build
                   )