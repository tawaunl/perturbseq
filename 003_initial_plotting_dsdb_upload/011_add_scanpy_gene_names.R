## Forgot to add the scanpy-based gene names in the initial
## upload. Add them here for convenience

library(magrittr)

library(Matrix)

library(SummarizedExperiment)
library(SingleCellExperiment)
library(MultiAssayExperiment)

library(dsassembly)

## Load unfiltered MAE

options(timeout=10000)

mae <- getDataset("DS000018385", version=1)
sce <- mae[['cellbender']]

## Load filtered SCE

initial_processing_config <- yaml::read_yaml('../002_initial_processing/config.yaml')
sce_filt <- readRDS(file.path(initial_processing_config$scratch_dir, 'convert_mdata_to_sce_filt.Rds'))

## Deal with duplicate symbols

symbol_id <- paste0(rowData(sce)$symbol, '_', rowData(sce)$gene_id)

old_name <- dplyr::if_else(
  symbol_id %in% rownames(sce_filt),
  symbol_id,
  rowData(sce)$symbol
)

stopifnot(!duplicated(old_name))
stopifnot(setequal(old_name, rownames(sce_filt)))

rowData(sce)$old_name <- old_name

## Add the names to datasetdb

mae[['cellbender']] <- sce
updateDataset(mae)

## Version 2 of data uploaded
