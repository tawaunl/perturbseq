library(magrittr)
library(Matrix)

library(SummarizedExperiment)
library(SingleCellExperiment)
library(MultiAssayExperiment)

library(reticulate)
use_python(python=Sys.which("python"))

muon <- import("muon")

#in_mdata <- '../scratch/output_3/add_wnn_leiden.h5mu'

args <- commandArgs(trailingOnly=TRUE)
in_mdata <- args[1]
out_rds <- args[2]

mdata <- muon$read(in_mdata)

# convert count matrices bcuz reticulate can't handle sparse int
# matrices
py_run_string("
X_gex_counts = r.mdata['GEX'].raw.X.astype(float)
X_gex_uncleaned_counts = r.mdata['GEX'].layers['uncleaned_counts'].astype(float)
X_adt_counts = r.mdata['ADT'].raw.X.astype(float)
X_hto_counts = r.mdata['HTO'].raw.X.astype(float)
X_sgrna_indicator = r.mdata['sgRNA'].X.astype(float)
X_sgrna_counts = r.mdata['sgRNA'].raw.X.astype(float)

has_wnn = 'connectivities' in r.mdata.obsp
has_wnn_leiden = 'leiden' in r.mdata.obs
has_wnn_umap = 'X_umap' in r.mdata.obsm
has_gex_nn = 'connectivities' in r.mdata['GEX'].obsp
has_adt_nn = 'connectivities' in r.mdata['ADT'].obsp
")

X_gex_logcounts <- as(t(mdata[['GEX']]$X), 'CsparseMatrix')
colnames(X_gex_logcounts) <- rownames(mdata[['GEX']]$obs)
rownames(X_gex_logcounts) <- rownames(mdata[['GEX']]$var)

X_gex_counts <- as(t(py$X_gex_counts), 'CsparseMatrix')
colnames(X_gex_counts) <- rownames(mdata[['GEX']]$obs)
rownames(X_gex_counts) <- rownames(mdata[['GEX']]$var)

X_gex_uncleaned_counts <- as(t(py$X_gex_uncleaned_counts), 'CsparseMatrix')
colnames(X_gex_uncleaned_counts) <- rownames(mdata[['GEX']]$obs)
rownames(X_gex_uncleaned_counts) <- rownames(mdata[['GEX']]$var)

X_adt_logcounts <- as(t(mdata[['ADT']]$X), 'CsparseMatrix')
colnames(X_adt_logcounts) <- rownames(mdata[['ADT']]$obs)
rownames(X_adt_logcounts) <- rownames(mdata[['ADT']]$var)

X_adt_counts <- as(t(py$X_adt_counts), 'CsparseMatrix')
colnames(X_adt_counts) <- rownames(mdata[['ADT']]$obs)
rownames(X_adt_counts) <- rownames(mdata[['ADT']]$var)

X_hto_counts <- as(t(py$X_hto_counts), 'CsparseMatrix')
colnames(X_hto_counts) <- rownames(mdata[['HTO']]$obs)
rownames(X_hto_counts) <- rownames(mdata[['HTO']]$var)

X_sgrna_counts <- as(t(py$X_sgrna_counts), 'CsparseMatrix')
colnames(X_sgrna_counts) <- rownames(mdata[['sgRNA']]$obs)
rownames(X_sgrna_counts) <- rownames(mdata[['sgRNA']]$var)

X_sgrna_indicator <- as(t(py$X_sgrna_indicator), 'CsparseMatrix')
colnames(X_sgrna_indicator) <- rownames(mdata[['sgRNA']]$obs)
rownames(X_sgrna_indicator) <- rownames(mdata[['sgRNA']]$var)

sce <- SingleCellExperiment(
  assays=list(logcounts=X_gex_logcounts,
              counts=X_gex_counts,
              rawcounts=X_gex_uncleaned_counts),
  rowData=mdata[['GEX']]$var,
  colData=mdata[['GEX']]$obs
)

# FIXME: Below assumes that mdata is in same order as
# mdata[GEX]. Technically it might be possible for mdata.obs to be in
# different order than mdata[GEX].obs, but shouldn't happen in
# practice. This assertion checks for that.
stopifnot(rownames(mdata$obs) == colnames(sce))

if (py$has_wnn_leiden) {
  colData(sce)$leiden_wnn <- mdata$obs$leiden
  colData(sce)$leiden_wnn2 <- mdata$obs$leiden2
}

if (py$has_wnn) {
  X_wnn <- as(t(mdata$obsp[['connectivities']]), 'CsparseMatrix')
  #rownames(X_wnn) <- colnames(X_wnn) <- rownames(mdata$obs)
  colPairs(sce)[['gex_adt_wnn_graph']] <- X_wnn
}

if (py$has_gex_nn) {
  X_gex_nn <- as(t(mdata[['GEX']]$obsp[['connectivities']]), 'CsparseMatrix')
  #rownames(X_gex_nn) <- colnames(X_gex_nn) <- rownames(mdata[['GEX']]$obs)
  colPairs(sce)[['gex_only_nn_graph']] <- X_gex_nn
}

get_adata_pca <- function(adata) {
  X_pca <- adata$obsm[['X_pca']]
  rownames(X_pca) <- rownames(adata$obs)
  colnames(X_pca) <- paste0('PC', 1:ncol(X_pca))

  attr(X_pca, "varExplained") <- adata$uns[['pca']]$variance

  # CHECK: Is it correct to multiply by 100?
  attr(X_pca, "percentVar") <- 100*adata$uns[['pca']]$variance_ratio

  # CHECK: Should I be using the transpose of this?
  attr(X_pca, "rotation") <- adata$varm$get('PCs')
  colnames(attr(X_pca, "rotation")) <- paste0('PC', 1:ncol(X_pca))
  rownames(attr(X_pca, "rotation")) <- rownames(adata$var)

  X_pca
}

get_adata_umap <- function(adata) {
  X_umap <- adata$obsm[['X_umap']]
  rownames(X_umap) <- rownames(adata$obs)
  colnames(X_umap) <- paste0('UMAP', 1:2)

  X_umap
}

reducedDim(sce, 'X_pca') <- get_adata_pca(mdata[['GEX']])
reducedDim(sce, 'X_umap') <- get_adata_umap(mdata[['GEX']])

if (py$has_wnn_umap) {
  reducedDim(sce, 'X_wnn_umap') <- get_adata_umap(mdata)
}

sce_hto <- SingleCellExperiment(
  assays=list(counts=X_hto_counts),
  rowData=mdata[['HTO']]$var,
  colData=mdata[['HTO']]$obs
)

reducedDim(sce_hto, 'X_umap') <- get_adata_umap(mdata[['HTO']])

sce_adt <- SingleCellExperiment(
  assays=list(logcounts=X_adt_logcounts, counts=X_adt_counts),
  rowData=mdata[['ADT']]$var,
  colData=mdata[['ADT']]$obs
)

reducedDim(sce_adt, 'X_pca') <- get_adata_pca(mdata[['ADT']])
reducedDim(sce_adt, 'X_umap') <- get_adata_umap(mdata[['ADT']])

if (py$has_adt_nn) {
  X_adt_nn <- as(t(mdata[['ADT']]$obsp[['connectivities']]), 'CsparseMatrix')
  #rownames(X_adt_nn) <- colnames(X_adt_nn) <- rownames(mdata[['GEX']]$obs)
  colPairs(sce_adt)[['nn_graph']] <- X_adt_nn
}

sce_sgrna <- SingleCellExperiment(
  assays=list(indicator=X_sgrna_indicator,
              counts=X_sgrna_counts),
  rowData=mdata[['sgRNA']]$var,
  colData=mdata[['sgRNA']]$obs
)

altExp(sce, 'Antibody Capture') <- sce_adt
altExp(sce, 'CRISPR Guide Capture') <- sce_sgrna
altExp(sce, 'Multiplexing Capture') <- sce_hto

mainExpName(sce) <- 'cellbender gene expression'

saveRDS(
  sce, out_rds,
  compress = FALSE
)
