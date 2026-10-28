# STEP 3: limma-trend differential expression, hypoxia vs normoxia
log <- file(snakemake@log[[1]], open = "wt"); sink(log); sink(log, type = "message")
suppressPackageStartupMessages(library(limma))

e  <- as.matrix(read.csv(snakemake@input[["expr"]], row.names = 1, check.names = FALSE))
sd <- read.csv(snakemake@input[["samples"]])
condition <- factor(sd$condition[match(colnames(e), sd$sample)], levels = c("normoxia", "hypoxia"))

fit <- eBayes(lmFit(e, model.matrix(~ condition)), trend = TRUE)
res <- topTable(fit, coef = "conditionhypoxia", number = Inf)
res <- data.frame(gene = rownames(res), log2FC = res$logFC, AveLog2 = res$AveExpr,
                  pvalue = res$P.Value, padj = res$adj.P.Val)
padj_cut <- snakemake@params[["padj"]]; lfc_cut <- snakemake@params[["lfc"]]
sig <- res[!is.na(res$padj) & res$padj < padj_cut & abs(res$log2FC) > lfc_cut, ]
write.csv(res, snakemake@output[["res"]], row.names = FALSE)
write.csv(sig, snakemake@output[["sig"]], row.names = FALSE)

hif <- c("CA9","VEGFA","SLC2A1","PDK1","BNIP3","NDRG1","EGLN3","P4HA1","ADM","ANKRD37",
         "LDHA","PGK1","BHLHE40","ANGPTL4","PFKFB3","HK2","DDIT4","NDUFA4L2","KDM3A","FAM162A")
up_hif <- intersect(hif, sig$gene[sig$log2FC > 0])
txt <- c("LNCaP hypoxia vs normoxia (GEO GSE106305, empty-vector cells)",
         sprintf("Genes tested: %d", nrow(res)),
         sprintf("Significant (padj < %s, |log2FC| > %s): %d  (up %d, down %d)", padj_cut, lfc_cut,
                 nrow(sig), sum(sig$log2FC > 0), sum(sig$log2FC < 0)),
         sprintf("Known HIF targets significantly up: %d of %d: %s", length(up_hif), length(hif),
                 paste(up_hif, collapse = ", ")),
         "", "Top 15 genes:", capture.output(print(head(res[, c("gene", "log2FC", "padj")], 15), row.names = FALSE)))
writeLines(txt, snakemake@output[["summary"]])
