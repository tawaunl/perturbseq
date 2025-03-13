library(magrittr)

library(Matrix)

library(SummarizedExperiment)
library(SingleCellExperiment)
library(MultiAssayExperiment)

library(dsassembly)

ANALYSIS_SCRATCH <- "/gstore/scratch/u/lucast3/ngs5388_collaborative/scratch"
mae <- getDataset("DS000018385")
sce <- mae[['cellbender']]
rownames(sce) <- make.unique(rowData(sce)$symbol)

# load assays to memory

convert_mat <- function(x) as(x, 'CsparseMatrix')
assays(sce) %<>% lapply(convert_mat)
assays(altExp(sce, 'Antibody Capture')) %<>% lapply(convert_mat)
assays(altExp(sce, 'CRISPR Guide Capture')) %<>% lapply(convert_mat)
assays(altExp(sce, 'Multiplexing Capture')) %<>% lapply(convert_mat)

# remove metadata (can cause issues in non-cedar environments, e.g. conda)

metadata(sce) <- list()
metadata(altExp(sce, 'Antibody Capture')) <- list()
metadata(altExp(sce, 'CRISPR Guide Capture')) <- list()
metadata(altExp(sce, 'Multiplexing Capture')) <- list()

# save it

dir.create(ANALYSIS_SCRATCH, recursive = TRUE)

saveRDS(
  sce,
  file.path(ANALYSIS_SCRATCH, 'converted_sce.Rds'),
  compress = FALSE
)
