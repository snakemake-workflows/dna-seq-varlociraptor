library(siglasso)
library(tibble)

# Load COSMIC signatures
cosmic_signatures <- read.table(snakemake@input[[1]], header = TRUE, stringsAsFactors = FALSE)
row.names(cosmic_signatures) <- cosmic_signatures[, 1]
cosmic_signatures <- cosmic_signatures[, -1]
row.names(cosmic_signatures) <- gsub("\\[|\\]", "", row.names(cosmic_signatures))
cosmic_signatures <- as.matrix(cosmic_signatures)

# Load sample context
first_line <- readLines(snakemake@input[[2]], n = 1)

if (length(first_line) != 0) {
    context_file <- read.table(snakemake@input[[2]])
    spectrum <- context2spec(context_file, plot=FALSE)

    sample_signatures <- rownames_to_column(
        as.data.frame(
            siglasso(spectrum, cosmic_signatures, plot=FALSE)
            ),
            var="Signature"
        )
} else {
    sample_signatures <- data.frame(matrix(ncol=2, nrow=0))
    colnames(sample_signatures) <- c("Signature", snakemake@wildcards[["group"]])
}
write.table(sample_signatures, file=snakemake@output[[1]], quote=FALSE, sep='\t')