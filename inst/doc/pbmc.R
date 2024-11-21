## ----include = FALSE----------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

## ----load-packages, results='hide', message=FALSE-----------------------------
# Load required packages
library(Seurat)
library(scregclust)

## ----download-data------------------------------------------------------------
url <- paste0(
  "https://cf.10xgenomics.com/samples/cell-arc/2.0.0/",
  "pbmc_granulocyte_sorted_3k/",
  "pbmc_granulocyte_sorted_3k_filtered_feature_bc_matrix.h5"
)
data_path <- file.path(
  tempdir(), "pbmc_granulocyte_sorted_3k_filtered_feature_bc_matrix.h5"
)

download.file(url, data_path, cacheOK = FALSE, mode = "wb")

## ----load-h5------------------------------------------------------------------
pbmc_data <- Read10X_h5(
  data_path,
  use.names = TRUE,
  unique.features = TRUE
)[["Gene Expression"]]

## ----create-seurat-object-----------------------------------------------------
pbmc <- CreateSeuratObject(
  counts = pbmc_data, min.cells = 3, min.features = 200
)

pbmc[["percent.mt"]] <- PercentageFeatureSet(pbmc, pattern = "^MT.")
pbmc <- subset(pbmc, subset = percent.mt < 30 & nFeature_RNA < 6000)

## ----apply-var-stabilization--------------------------------------------------
pbmc <- SCTransform(pbmc, variable.features.n = 6000)

z <- GetAssayData(pbmc, layer = "scale.data")
dim(z)

## ----prep-scregclust----------------------------------------------------------
out <- scregclust_format(z, mode = "TF")

## ----extract-scregclust-arguments---------------------------------------------
genesymbols <- out$genesymbols
sample_assignment <- out$sample_assignment
is_regulator <- out$is_regulator

## ----run-scregclust-----------------------------------------------------------
# set.seed(8374)
# fit <- scregclust(
#   z, genesymbols, is_regulator, penalization = seq(0.1, 0.5, 0.05),
#   n_modules = 10L, n_cycles = 50L, noise_threshold = 0.05
# )
# saveRDS(fit, file = "pbmc_scregclust.rds")

url <- paste0(
  "https://github.com/scmethods/scregclust/raw/main/datasets/",
  "pbmc_scregclust.rds"
)
fit_path <- file.path(tempdir(), "pbmc_scregclust.rds")
download.file(url, fit_path)
fit <- readRDS(fit_path)

## ----viz-metrics, fig.width=7, fig.height=4, fig.dpi=100----------------------
plot(fit)

## ----n-configs----------------------------------------------------------------
sapply(fit$results, function(r) length(r$output))

## ----viz-reg-network, fig.width=7, fig.height=7, fig.dpi=100------------------
plot_regulator_network(fit$results[[1]]$output[[1]])

