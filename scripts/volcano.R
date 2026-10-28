# STEP 4b: volcano plot
log <- file(snakemake@log[[1]], open = "wt"); sink(log); sink(log, type = "message")
suppressPackageStartupMessages(library(ggplot2))
res <- read.csv(snakemake@input[[1]])
padj <- snakemake@params[["padj"]]; lfc <- snakemake@params[["lfc"]]
res$status <- "NS"
res$status[res$padj < padj & res$log2FC >  lfc] <- "Up"
res$status[res$padj < padj & res$log2FC < -lfc] <- "Down"
hif <- c("CA9","VEGFA","SLC2A1","PDK1","BNIP3","NDRG1","EGLN3","P4HA1","ADM","ANKRD37")
lab <- res[res$gene %in% hif & res$status == "Up", ]
if (nrow(lab) < 3) lab <- head(res[res$status != "NS", ], 10)
g <- ggplot(res, aes(log2FC, -log10(padj), colour = status)) + geom_point(size = .8, alpha = .6) +
  scale_colour_manual(values = c(Down = "#2166ac", NS = "grey78", Up = "#b2182b")) +
  geom_vline(xintercept = c(-lfc, lfc), linetype = "dashed", colour = "grey40") +
  geom_hline(yintercept = -log10(padj), linetype = "dashed", colour = "grey40") +
  geom_text(data = lab, aes(label = gene), colour = "black", size = 3.3, vjust = -.8, fontface = "bold") +
  labs(title = "Volcano: hypoxia vs normoxia", subtitle = "GEO GSE106305", x = "log2 fold change",
       y = "-log10 adjusted p", colour = NULL) + theme_bw(base_size = 13) + theme(legend.position = "bottom")
ggsave(snakemake@output[[1]], g, width = 6.5, height = 6, dpi = 200)
