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
reticulate::use_python("~/.conda/envs/MOFA/bin/python", required = T)
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

Perturbations=as.matrix(as(assay(sce_crispr, 'counts'),
                           'CsparseMatrix'))


# indicator matrix pooling all guides within a gene KO
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


# Create MOFA object ---------


MOFAobject <- create_mofa(data)

# Visualise data structure
cairo_pdf("DataOverviewGuideLevel.pdf")
plot_data_overview(MOFAobject)
dev.off()

# Define options----------------------


## Data options ----------
# - scale_views: if views have very different ranges/variances, it is good practice to scale each view to unit variance (default is FALSE)
data_opts <- get_default_data_options(MOFAobject)
data_opts$scale_views <- TRUE
## Model options -------------
## - likelihoods: likelihood per view (options are "gaussian","poisson","bernoulli"). "gaussian" is used by default
## - num_factors: number of factors. By default K=10
model_opts <- get_default_model_options(MOFAobject)
model_opts$num_factors <- 12

## Training options----=---------
## - maxiter: number of iterations
## - convergence_mode: "fast", "medium", "slow". For exploration, the fast mode is good enough.
## - drop_factor_threshold: minimum variance explained criteria to drop factors while training. Default is -1 (no dropping of factors)
## - gpu_mode: use GPU mode? This needs cupy installed and a functional GPU, see https://biofam.github.io/MOFA2/gpu_training.html
## - seed: random seed
train_opts <- get_default_training_options(MOFAobject)
train_opts$convergence_mode <- "fast"
train_opts$seed <- 824
train_opts$drop_factor_threshold <- -1
#Prepare MOFA object ----------------

MOFAobject <- prepare_mofa(MOFAobject,
                           data_options = data_opts,
                           training_options = train_opts,
                           model_options = model_opts)


# Train the model-----------------------------


MOFAobject <- run_mofa(MOFAobject)

# Save the model --------------

outfile <- paste0(getwd(),"/model_GuideLevelCountswCorrectADT.hdf5")
saveRDS(MOFAobject, outfile)
