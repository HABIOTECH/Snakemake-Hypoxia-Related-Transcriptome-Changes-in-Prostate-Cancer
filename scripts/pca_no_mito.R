# Sensitivity analysis: PCA excluding mitochondrial genes
log <- file(snakemake@log[[1]], open = "wt")
sink(log)
sink(log, type = "message")

suppressPackageStartupMessages(library(ggplot2))

e <- as.matrix(read.csv(
  snakemake@input[["expr"]],
  row.names = 1,
  check.names = FALSE
))

sd <- read.csv(snakemake@input[["samples"]])

# Exclude mitochondrial genes by HGNC-style gene symbol
is_mito <- grepl("^MT-", rownames(e), ignore.case = TRUE)
cat("Genes before mitochondrial exclusion:", nrow(e), "\n")
cat("Mitochondrial genes excluded:", sum(is_mito), "\n")

e <- e[!is_mito, , drop = FALSE]
stopifnot(nrow(e) >= 2)

# Match the original PCA: select the 500 most variable genes
v <- apply(e, 1, var)
e <- e[order(-v)[seq_len(min(500, nrow(e)))], , drop = FALSE]

pc <- prcomp(t(e))
pv <- round(100 * pc$sdev^2 / sum(pc$sdev^2))

pd <- data.frame(
  PC1 = pc$x[, 1],
  PC2 = pc$x[, 2],
  sample = colnames(e),
  condition = factor(
    sd$condition[match(colnames(e), sd$sample)],
    levels = c("normoxia", "hypoxia")
  )
)

p <- ggplot(
  pd,
  aes(PC1, PC2, colour = condition, label = sample)
) +
  geom_point(size = 4) +
  geom_text(vjust = -1.3, show.legend = FALSE) +
  scale_colour_manual(
    values = c(normoxia = "#2166ac", hypoxia = "#b2182b")
  ) +
  scale_x_continuous(expand = expansion(mult = .25)) +
  scale_y_continuous(expand = expansion(mult = .3)) +
  labs(
    title = "PCA excluding mitochondrial genes",
    subtitle = "LNCaP hypoxia vs normoxia | GEO GSE106305",
    x = paste0("PC1 (", pv[1], "%)"),
    y = paste0("PC2 (", pv[2], "%)"),
    colour = NULL
  ) +
  theme_bw(base_size = 13) +
  theme(legend.position = "bottom")

ggsave(
  snakemake@output[[1]],
  p,
  width = 6,
  height = 5.5,
  dpi = 200
)

sink(type = "message")
sink()
close(log)
