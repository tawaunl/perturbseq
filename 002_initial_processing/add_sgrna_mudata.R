library(magrittr)
library(Matrix)
library(SummarizedExperiment)

library(reticulate)

use_python(python=Sys.which("python"))

anndata <- import("anndata")
muon <- import("muon")

#in_mdata <- "scratch/231030_apply_qc_fix_trem2.h5mu"
#in_se <- "scratch/231030_make_sgrna_summarized_experiment.Rds"
#out_mdata <- "scratch/231030_add_sgrna_mudata.h5mu"
args <- commandArgs(trailingOnly=TRUE)
in_mdata <- args[1]
in_se <- args[2]
out_mdata <- args[3]

sgrna_se <- readRDS(in_se)
mdata <- muon$read(in_mdata)

# Subset the sgRNA UMI count matrix to the valid cells

mat_sub <- assays(sgrna_se)[['union']][rownames(mdata[['GEX']]$obs),]

# get some cell-level statistics

data.frame(
  n_umi = rowSums(mat_sub)
) %>%
  dplyr::mutate(simpson_divers=(n_umi^2)/rowSums(mat_sub^2)) %>%
  dplyr::mutate(
    max_umi_cnt=sparseMatrixStats::rowMaxs(mat_sub),
    second_umi_cnt=sparseMatrixStats::rowOrderStats(mat_sub, which=ncol(mat_sub)-1)
  ) %>%
  dplyr::mutate(
    max_umi_frac = max_umi_cnt / n_umi
  ) ->
  meta

# get the top guide for each cell (NA if there's a tie)

mat_sub_Tsparse <- as(mat_sub, 'TsparseMatrix')

data.frame(
  i=mat_sub_Tsparse@i+1,
  j=mat_sub_Tsparse@j+1,
  x=mat_sub_Tsparse@x
) %>%
  dplyr::filter(
    meta$max_umi_cnt[i] > meta$second_umi_cnt[i],
    x == meta$max_umi_cnt[i]
  ) ->
  df

stopifnot(!duplicated(df$i))

top_guide <- rep(NA_character_, nrow(meta))
top_guide[df$i] <- colnames(mat_sub)[df$j]

meta$top_sgrna <- top_guide

# sanity checks
x <- mat_sub[cbind(
  1:nrow(mat_sub),
  as.integer(factor(top_guide,
                    levels=colnames(mat_sub)))
)[!is.na(top_guide),]]

stopifnot(x > 0)
stopifnot(x == sparseMatrixStats::rowMaxs(mat_sub)[!is.na(top_guide)])

# make the sgRNA anndata

sgrna_var <- data.frame(
  total_umis=colSums(mat_sub)
)
rownames(sgrna_var) <- colnames(mat_sub)

mat_sub_Rsparse <- as(mat_sub, 'RsparseMatrix')

adata = anndata$AnnData(
  X=mat_sub_Rsparse,
  obs=meta,
  var=sgrna_var
)

# make the new mudata and save it

mdata = muon$MuData(list(
  GEX=mdata[['GEX']],
  ADT=mdata[['ADT']],
  HTO=mdata[['HTO']],
  sgRNA=adata
))

mdata$write(out_mdata)
