# Sensitivity analysis: limma-trend excluding mitochondrial genes
suppressPackageStartupMessages(library(limma))

e <- as.matrix(read.csv(
    snakemake@input[["expr"]],
    row.names = 1,
    check.names = FALSE
))

sd <- read.csv(snakemake@input[["samples"]])

# Exclude mitochondrial genes (human gene symbols beginning MT-)
keep <- !grepl("^MT-", rownames(e), ignore.case = TRUE)
e <- e[keep, , drop = FALSE]

condition <- factor(
    sd$condition[match(colnames(e), sd$sample)],
    levels = c("normoxia", "hypoxia")
)

fit <- eBayes(lmFit(e, model.matrix(~ condition)), trend = TRUE)
res <- topTable(fit, coef = "conditionhypoxia", number = Inf)

res <- data.frame(
    gene = rownames(res),
    log2FC = res$logFC,
    AveLog2 = res$AveExpr,
    pvalue = res$P.Value,
    padj = res$adj.P.Val
)

padj_cut <- snakemake@params[["padj"]]
lfc_cut <- snakemake@params[["lfc"]]

sig <- res[
    !is.na(res$padj) &
    res$padj < padj_cut &
    abs(res$log2FC) > lfc_cut,
]

write.csv(res, snakemake@output[["res"]], row.names = FALSE)
write.csv(sig, snakemake@output[["sig"]], row.names = FALSE)

cat("Mitochondrial genes excluded:", sum(!keep), "\n")
cat("Genes tested:", nrow(res), "\n")
cat("Significant genes:", nrow(sig), "\n")
cat("Upregulated:", sum(sig$log2FC > 0), "\n")
cat("Downregulated:", sum(sig$log2FC < 0), "\n")
