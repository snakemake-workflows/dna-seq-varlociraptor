library(siglasso)
library(tibble)
library(dplyr)
library(readr)


# Load COSMIC signatures
cosmic_signatures <- read_tsv(snakemake@input[[1]])
cosmic_signatures <-  as.matrix(cosmic_signatures
                                %>% mutate(Type = gsub("\\[|\\]", "", Type))
                                %>% column_to_rownames(var = "Type")
                    )

sample_substitutions <- read_tsv(snakemake@input[[2]])
if (nrow(sample_substitutions) == 0) {
    for (output_file in snakemake@output) {
        write_tsv(tibble(), output_file)
    }
} else {
    for (output_file in snakemake@output) {
        min_vaf <- as.numeric(strsplit(output_file, split="\\.")[[1]][3]) / 100
        filtered_substitions <- (
            sample_substitutions
            %>% filter(AF >= min_vaf)
            %>% mutate(AF = NULL)
        )
        if (nrow(filtered_substitions) == 0) {
            write_tsv(tibble(), output_file)
        } else {
            spectrum <- context2spec(filtered_substitions, plot=FALSE)
            sample_signatures <- (
                as.data.frame(siglasso(spectrum, cosmic_signatures, plot=FALSE)) 
                %>% rownames_to_column(var="Signature")
                %>% filter(!!sym(snakemake@wildcards[["group"]]) > 0)
                %>% add_column(Frequency = min_vaf)
            )
            write_tsv(sample_signatures, output_file, col_names=FALSE)
        }
    }
}