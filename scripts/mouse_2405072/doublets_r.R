#!/usr/bin/env Rscript
# scDblFinder + DoubletFinder for one sample (called by l1_mouse.py --stage doublets).
# Usage: Rscript doublets_r.R <10x_dir> <barcodes.txt> <expected_rate> <out.csv>
# Each method is optional: if its package is missing or it fails, its columns are NA and the
# Python side drops it from the consensus (and logs that).
suppressPackageStartupMessages(library(Matrix))
args <- commandArgs(trailingOnly = TRUE)
dir <- args[1]; bc_file <- args[2]; rate <- as.numeric(args[3]); out <- args[4]

f <- function(stem) { p <- file.path(dir, stem); if (file.exists(paste0(p, ".gz"))) paste0(p, ".gz") else p }
m <- readMM(f("matrix.mtx"))
feat <- read.delim(f("features.tsv"), header = FALSE, stringsAsFactors = FALSE)
bcs <- read.delim(f("barcodes.tsv"), header = FALSE, stringsAsFactors = FALSE)[, 1]
if (ncol(feat) >= 3) {
  gex <- feat[, 3] == "Gene Expression"
  m <- m[gex, ]
  feat <- feat[gex, ]
}
rownames(m) <- make.unique(feat[, 2]); colnames(m) <- bcs
keep <- readLines(bc_file)
m <- as(m[, keep], "CsparseMatrix")
cat("cells", ncol(m), "genes", nrow(m), "rate", rate, "\n")

res <- data.frame(barcode = keep, scdbl_score = NA_real_, scdbl_call = NA, df_pann = NA_real_, df_call = NA)

# --- scDblFinder
tryCatch({
  suppressPackageStartupMessages({ library(SingleCellExperiment); library(scDblFinder) })
  set.seed(1)
  sce <- scDblFinder(SingleCellExperiment(list(counts = m)), dbr = rate)
  res$scdbl_score <- sce$scDblFinder.score
  res$scdbl_call <- sce$scDblFinder.class == "doublet"
  cat("scDblFinder calls", sum(res$scdbl_call), "\n")
}, error = function(e) cat("scDblFinder failed:", conditionMessage(e), "\n"))

# --- DoubletFinder (handles the _v3 / new function names)
tryCatch({
  suppressPackageStartupMessages({ library(Seurat); library(DoubletFinder) })
  set.seed(1)
  so <- CreateSeuratObject(m)
  so <- NormalizeData(so, verbose = FALSE)
  so <- FindVariableFeatures(so, nfeatures = 2000, verbose = FALSE)
  so <- ScaleData(so, verbose = FALSE)
  so <- RunPCA(so, npcs = 30, verbose = FALSE)
  so <- FindNeighbors(so, dims = 1:20, verbose = FALSE)
  so <- FindClusters(so, resolution = 0.5, verbose = FALSE)
  sweep_fn <- if (exists("paramSweep")) paramSweep else paramSweep_v3
  df_fn <- if (exists("doubletFinder")) doubletFinder else doubletFinder_v3
  sw <- summarizeSweep(sweep_fn(so, PCs = 1:20, sct = FALSE), GT = FALSE)
  bc <- find.pK(sw)
  pK <- as.numeric(as.character(bc$pK[which.max(bc$BCmetric)]))
  homo <- modelHomotypic(so$seurat_clusters)
  nexp <- round(rate * ncol(so) * (1 - homo))
  so <- tryCatch(df_fn(so, PCs = 1:20, pN = 0.25, pK = pK, nExp = nexp, reuse.pANN = NULL, sct = FALSE),
                 error = function(e) df_fn(so, PCs = 1:20, pN = 0.25, pK = pK, nExp = nexp, reuse.pANN = FALSE, sct = FALSE))
  md <- so@meta.data
  pc <- grep("^pANN", colnames(md), value = TRUE)[1]
  dc <- grep("^DF.classifications", colnames(md), value = TRUE)[1]
  res$df_pann <- md[keep, pc]
  res$df_call <- md[keep, dc] == "Doublet"
  cat("DoubletFinder pK", pK, "calls", sum(res$df_call), "\n")
}, error = function(e) cat("DoubletFinder failed:", conditionMessage(e), "\n"))

write.csv(res, out, row.names = FALSE)
