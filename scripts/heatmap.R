# STEP 4c: heatmap of the top 30 genes
log <- file(snakemake@log[[1]], open = "wt"); sink(log); sink(log, type = "message")
suppressPackageStartupMessages(library(pheatmap))
e   <- as.matrix(read.csv(snakemake@input[["expr"]], row.names = 1, check.names = FALSE))
sd  <- read.csv(snakemake@input[["samples"]])
res <- read.csv(snakemake@input[["res"]])
top <- head(res$gene[res$gene %in% rownames(e)], 30)
m <- e[top, ]; m <- m - rowMeans(m)
ann <- data.frame(condition = sd$condition[match(colnames(m), sd$sample)], row.names = colnames(m))
pheatmap(m, annotation_col = ann,
         annotation_colors = list(condition = c(normoxia = "#2166ac", hypoxia = "#b2182b")),
         main = "Top 30 genes (GSE106305)", border_color = NA, angle_col = 0,
         filename = snakemake@output[[1]], width = 6, height = 8)
