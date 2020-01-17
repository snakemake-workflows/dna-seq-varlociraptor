library(ComplexHeatmap)
library(ggplot2)

table = read.table(snakemake@input[[1]], sep="\t", header=TRUE, row.names=1)
mat = as.matrix(table)
mat = t(mat)


col = c(SNV = "blue", INDEL = "red")

alter_fun = list(
        SNV = function(x, y, w, h) grid.rect(x, y, w*0.9, h*0.9, 
            gp = gpar(fill = col["SNV"], col = NA)),
        INDEL = function(x, y, w, h) grid.rect(x, y, w*0.9, h*0.4, 
            gp = gpar(fill = col["INDEL"], col = NA))
    )


heatmap_legend_param = list(title = "Alterations", at = c("SNV", "INDEL"), 
        labels = c("SNV", "INDEL"))

mat <- mat[order(apply(mat, 1, function(row) sum(row != "")), decreasing = T), ]

if (nrow(mat) > 2000) {
    mat <- mat[1:2000,]
}
rows_matrix <- nrow(mat)
height_plot <- (rows_matrix/5)
if (height_plot < 4) {
    height_plot <- 4
}
pdf(file = snakemake@output[[1]], height=height_plot)
if (rows_matrix > 0) {
    oncoprint <- oncoPrint(mat,
        alter_fun = alter_fun, col = col, 
        remove_empty_columns = FALSE, remove_empty_rows = TRUE,
        pct_side = "right", row_names_side = "left",
        show_column_names=T,
        column_title = "OncoPrint", heatmap_legend_param = heatmap_legend_param)
    draw(oncoprint, newpage=F)
}
dev.off()