# =====================================================================
#  Snakemake pipeline: LNCaP prostate cancer cells, hypoxia vs normoxia
#  Real data: GEO GSE106305 (Guo et al., 2019, Nature Communications)
#
#  Run one step at a time:
#     snakemake --use-conda --cores 2 step1_download
#     snakemake --use-conda --cores 2 step2_prepare
#     snakemake --use-conda --cores 2 step3_de
#     snakemake --use-conda --cores 2 step4_plots
#  Or everything:
#     snakemake --use-conda --cores 2
# =====================================================================
configfile: "config/config.yaml"

# ----------------------------- targets -------------------------------
rule all:
    input:
        "results/summary.txt",
        "results/PCA.png",
        "results/volcano.png",
        "results/heatmap_top30.png"

rule step1_download:
    input: "data/GSE106305_LNCaP_RPKM.txt.gz"

rule step2_prepare:
    input: "results/expression_log2.csv"

rule step3_de:
    input: "results/DE_results.csv", "results/summary.txt"

rule step4_plots:
    input: "results/PCA.png", "results/volcano.png", "results/heatmap_top30.png"

# ----------------------- STEP 1: download data ------------------------
rule download:
    output: "data/GSE106305_LNCaP_RPKM.txt.gz"
    params: url = config["geo_url"]
    conda: "envs/download.yaml"
    log: "logs/download.log"
    shell: "curl -fsSL -o {output} {params.url} 2> {log}"

# ------------- STEP 2: pick samples, filter, log-transform -------------
rule prepare:
    input: "data/GSE106305_LNCaP_RPKM.txt.gz"
    output:
        expr    = "results/expression_log2.csv",
        samples = "results/samples.csv"
    params:
        samples  = config["samples"],
        min_rpkm = config["min_rpkm"]
    conda: "envs/r.yaml"
    log: "logs/prepare.log"
    script: "scripts/prepare.R"

# ------------- STEP 3: differential expression (limma) ----------------
rule de:
    input:
        expr    = "results/expression_log2.csv",
        samples = "results/samples.csv"
    output:
        res     = "results/DE_results.csv",
        sig     = "results/DE_significant.csv",
        summary = "results/summary.txt"
    params:
        padj = config["padj_cutoff"],
        lfc  = config["lfc_cutoff"]
    conda: "envs/r.yaml"
    log: "logs/de.log"
    script: "scripts/de.R"

# ------------------------- STEP 4: plots ------------------------------
rule pca:
    input:
        expr    = "results/expression_log2.csv",
        samples = "results/samples.csv"
    output: "results/PCA.png"
    conda: "envs/r.yaml"
    log: "logs/pca.log"
    script: "scripts/pca.R"

rule pca_no_mito:
    input:
        expr = "results/expression_log2.csv",
        samples = "results/samples.csv"
    output: "results/PCA_no_mito.png"
    conda: "envs/r.yaml"
    log: "logs/pca_no_mito.log"
    script: "scripts/pca_no_mito.R"

rule volcano:
    input: "results/DE_results.csv"
    output: "results/volcano.png"
    params:
        padj = config["padj_cutoff"],
        lfc  = config["lfc_cutoff"]
    conda: "envs/r.yaml"
    log: "logs/volcano.log"
    script: "scripts/volcano.R"

rule heatmap:
    input:
        expr    = "results/expression_log2.csv",
        samples = "results/samples.csv",
        res     = "results/DE_results.csv"
    output: "results/heatmap_top30.png"
    conda: "envs/r.yaml"
    log: "logs/heatmap.log"
    script: "scripts/heatmap.R"

# -------- Sensitivity analysis: differential expression without MT genes --------
rule de_no_mito:
    input:
        expr    = "results/expression_log2.csv",
        samples = "results/samples.csv"
    output:
        res = "results/DE_results_no_mito.csv",
        sig = "results/DE_significant_no_mito.csv"
    params:
        padj = config["padj_cutoff"],
        lfc  = config["lfc_cutoff"]
    conda: "envs/r.yaml"
    log: "logs/de_no_mito.log"
    script: "scripts/de_no_mito.R"
