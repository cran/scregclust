## ----include = FALSE----------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

## ----load-packages, results='hide', message=FALSE-----------------------------
# Load required packages
library(GEOquery)
library(Seurat)
library(scregclust)

## ----load-data, results='hide', message=FALSE, warning=FALSE------------------
# Download the gene expression data
url <- paste0(
  "https://www.ncbi.nlm.nih.gov/geo/download/",
  "?acc=GSE60361&format=file&",
  "file=GSE60361%5FC1%2D3005%2DExpression%2Etxt%2Egz"
)
expr_path <- file.path(tempdir(), "Expression.txt.gz")
download.file(url, expr_path, cacheOK = FALSE, mode = "wb")

# Load the gene expression data
expr <- read.table(
  expr_path,
  header = TRUE,
  sep = "\t",
  stringsAsFactors = FALSE,
  fill = TRUE
)

# A few gene symbols appear as duplicates, make unique.
gene_symbols <- make.unique(expr[, 1], sep = "-")
expr <- expr[, -1]
rownames(expr) <- gene_symbols

# Download meta data
gse <- getGEO("GSE60361")
meta_data <- pData(phenoData(gse[[1]]))
# Sample names are stored in the meta data's row names
sample_names <- rownames(meta_data)
colnames(expr) <- sample_names

# Create Seurat object and preprocess the data using SCTransform
mouse <- CreateSeuratObject(
  counts = expr,
  min.cells = 3,
  min.features = 500,
  meta.data = meta_data
)
mouse <- SCTransform(mouse, verbose = TRUE)

## ----load-tfs, results='hide', message=FALSE----------------------------------
url <- "https://resources.aertslab.org/cistarget/tf_lists/allTFs_mm.txt"
tfs_path <- file.path(tempdir(), "allTFs_mm.txt")
download.file(url, tfs_path, cacheOK = FALSE, mode = "w")
tfs <- read.table(
  tfs_path,
  header = FALSE,
  sep = "\t",
  stringsAsFactors = FALSE
)
tfs <- tfs[, 1]

## ----extract-gene-cells-table-------------------------------------------------
z <- GetAssayData(mouse, layer = "scale.data")
dim(z)

## ----scregclust-format--------------------------------------------------------
out <- scregclust_format(z, mode = "TF")

genesymbols <- out$genesymbols
sample_assignment <- out$sample_assignment

## ----manual-is-regulator------------------------------------------------------
is_regulator <- rep(0, length = length(genesymbols))
is_regulator[which(genesymbols %in% tfs)] <- 1

## ----run-scregclust-----------------------------------------------------------
# # Run scregclust
# set.seed(8374)
# fit <- scregclust(
#   z, genesymbols, is_regulator, penalization = seq(0.1, 0.5, 0.05),
#   n_modules = 10L, n_cycles = 50L, noise_threshold = 0.05
# )
# saveRDS(fit, file = "datasets/mouse_scregclust.rds")

url <- paste0(
  "https://github.com/scmethods/scregclust/raw/main/datasets/",
  "mouse_scregclust.rds"
)
fit_path <- file.path(tempdir(), "mouse_scregclust.rds")
download.file(url, fit_path)
fit <- readRDS(fit_path)

## ----viz-fit, fig.width=7, fig.height=4, fig.dpi=100--------------------------
plot(fit)

