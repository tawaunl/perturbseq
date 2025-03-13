library(magrittr)

library(Matrix)

library(SummarizedExperiment)
library(SingleCellExperiment)
library(MultiAssayExperiment)

library(dsassembly)

## Load unfiltered MAE

options(timeout=10000)

# Fix issue with corrupted file in cache
# due to my session crashing on initial download
#gp.cache::updateCache(TRUE)

mae <- getDataset("DS000018385", version=2)
sce <- mae[['cellbender']]

## Load filtered SCE

initial_processing_config <- yaml::read_yaml('../002_initial_processing/config.yaml')
sce_filt <- readRDS(file.path(initial_processing_config$scratch_dir, 'convert_mdata_to_sce_filt.Rds'))

## Filter SCE, make the colnames align

stopifnot(colnames(sce_filt) %in% colnames(sce))
sce <- sce[, colnames(sce_filt)]

## Make rownames align, add GEX logcounts

stopifnot(!duplicated(rowData(sce)$old_name))
stopifnot(setequal(rowData(sce)$old_name, rownames(sce_filt)))

sce_filt <- sce_filt[rowData(sce)$old_name,]
rownames(sce_filt) <- rownames(sce)

logcounts(sce) <- logcounts(sce_filt)

## Update GEX reducedDims

reducedDims(sce) <- reducedDims(sce_filt)

## update ADT reducedDims

reducedDims(altExp(sce, 'Antibody Capture')) <- reducedDims(altExp(sce_filt, 'Antibody Capture'))

## Add leiden, leiden2

colData(sce)$leiden <- NULL
colData(altExp(sce, 'Antibody Capture'))$leiden <- NULL

colData(sce)$leiden1 <- colData(sce_filt)$leiden_wnn
colData(sce)$leiden2 <- colData(sce_filt)$leiden_wnn2

## Add sce back to mae

mae[['cellbender']] <- sce

## Update dataset

updateDataset(mae)

## Successfully uploaded version 3
