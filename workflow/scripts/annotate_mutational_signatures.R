library(siglasso)


# Load COSMIC signatures
cosmic_signatures <- read.table(snakemake@input[[1]], header = TRUE, stringsAsFactors = FALSE)
row.names(cosmic_signatures) <- cosmic_signatures[, 1]
cosmic_signatures <- cosmic_signatures[, -1]
row.names(cosmic_signatures) <- gsub("\\[|\\]", "", row.names(cosmic_signatures))
cosmic_signatures <- as.matrix(cosmic_signatures)

# Load sample context
context_file <- read.table(snakemake@input[[2]])
spectrum <- context2spec(context_file, plot=FALSE)

sample_signatures <- siglasso(spectrum, cosmic_signatures, plot=FALSE)
write.table(sample_signatures, file=snakemake@output[[1]], quote=FALSE, sep='\t')