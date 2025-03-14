library(MOFA2)
library(data.table)
library(magrittr)
library(dtplyr)
library(dplyr)
library(Matrix)
library(MultiAssayExperiment)
library(SingleCellExperiment)
library(dsassembly)
library(tidyverse)
# (Optional) set up reticulate connection with Python
# library(reticulate)
setwd("/gstore/scratch/u/lucast3/ngs5388_collaborative/MOFA")
########Load Data##############
## Load the SCE object and extract some objects
mae <- getDataset("DS000018385")
sce <- mae[['cellbender']]
rownames(sce) <- make.unique(rowData(sce)$symbol)

sce_crispr <- altExp(sce, 'CRISPR Guide Capture')
cite <- altExp(sce, 'Antibody Capture')

# grna by cell, 1=cell has grna, 0=else
grna_indicator <- as(assay(sce_crispr, 'indicator'),
                     'CsparseMatrix')

Perturbations=as.matrix(grna_indicator)

# indicator matrix pooling all guides within a gene KO

# Normalize data -----------
library(scater)
library(scran)
library(limma)
set.seed(824)  # For reproducibility
library(Seurat)
library(future)
options(future.globals.maxSize= 52428800000)

# Assuming your data is in a SingleCellExperiment (coded as 'sce'), you'll need to convert it
# Start by creating a Seurat object from SCE
counts(sce) <- as(counts(sce),"dgCMatrix")
sce_to_seurat <- CreateSeuratObject(counts = counts(sce),meta.data = data.frame(colData(sce)))

logcounts(cite) <- as(logcounts(cite),"dgCMatrix")

# Normalizing the data
plan("multicore", workers = 6)
sce_to_seurat <- NormalizeData(sce_to_seurat, normalization.method = "LogNormalize", scale.factor = 10000)
sce_to_seurat <- ScaleData(sce_to_seurat, vars.to.regress = c("n_genes_by_counts","total_counts"))

sce_to_seurat <- FindVariableFeatures(sce_to_seurat,nfeatures = 5000)



RNA <- GetAssayData(sce_to_seurat,slot = "scale.data")
RNA <- RNA[which(rownames(RNA) %in% VariableFeatures(sce_to_seurat)),]

ADT <- logcounts(cite)
vars <- data.frame(n_genes=cite$n_genes_by_counts,totalCounts=cite$total_counts)
ADT <- removeBatchEffect(ADT, covariates = vars, group = factor(sce$gate))

data <- list(RNA=as.matrix(RNA),
             Protein=as.matrix(ADT),
             Perturbations=Perturbations)
saveRDS(data,"ListofData.rds")

output_dir <- file.path(getwd(),"numpydata")

lapply(names(data), function(view_name) {
  write.csv(t(data[[view_name]]), file = file.path(output_dir, paste0(view_name, ".csv")), row.names = FALSE)
})
