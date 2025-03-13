library(magrittr)
library(Matrix)
library(SummarizedExperiment)

library(reticulate)
use_python(python=Sys.which("python"))
muon <- import("muon")

#in_barcodes <- "scratch/231009_merge_raw_multimodal_barcodes.txt.gz"
#in_counts <- "scratch/231007_combine_sgrna_counts.csv"
#in_mdata <- "scratch/231030_apply_qc_fix_trem2.h5mu"
#out_se <- "scratch/231030_make_sgrna_summarized_experiment.Rds"

args <- commandArgs(trailingOnly=TRUE)
in_barcodes <- args[1]
in_counts <- args[2]
in_mdata <- args[3]
out_se <- args[4]

mdata <- muon$read(in_mdata, backed=TRUE)

raw_cell_ids <- scan(
  in_barcodes,
  what='character'
)

stopifnot(rownames(mdata[['GEX']]$obs) %in% raw_cell_ids)

sgrna_counts_all <- read.csv(
  in_counts,
  stringsAsFactors = TRUE
)

sgrna_counts_all %<>%
  dplyr::mutate(
    cell_id=as.factor(paste0(SAMID, "_", Barcode))
  )

sgrna_counts_all %<>%
  dplyr::mutate(valid_cell=cell_id %in% rownames(mdata[['GEX']]$obs))

pass_qc <- !(mdata[['GEX']]$obs$is_outlier)
names(pass_qc) <- rownames(mdata[['GEX']]$obs)

sgrna_counts_all$pass_qc <- pass_qc[as.character(sgrna_counts_all$cell_id)]

n_htos_detected <- mdata[['GEX']]$obs$n_htos_detected
names(n_htos_detected) <- rownames(mdata[['GEX']]$obs)

sgrna_counts_all$n_htos <- n_htos_detected[as.character(sgrna_counts_all$cell_id)]

sgrna_counts_all %>%
  #dplyr::filter(valid_cell) %>%
  dplyr::select(SAMID, cell_id, portion, UMI, Feature, Count, Total,
                valid_cell, pass_qc, n_htos) %>%
  tidyr::pivot_wider(names_from=portion, values_from=c(Count, Total),
                     values_fill=0) %>%
  dplyr::mutate(Count_Both=Count_Pellet + Count_Supernatant,
                Total_Both=Total_Pellet + Total_Supernatant) ->
  sgrna_counts_pivot

sgrna_counts_pivot %>%
  dplyr::group_by(SAMID, valid_cell, pass_qc, n_htos, cell_id, Feature) %>%
  dplyr::summarize(pellet_only=sum(Count_Pellet > 0 & Count_Supernatant == 0),
                   supernatant_only=sum(Count_Pellet == 0 & Count_Supernatant > 0),
                   both=sum(Count_Pellet > 0 & Count_Supernatant > 0),
                   .groups='drop') ->
  sgrna_pivot_sum_umi

sgrna_pivot_sum_umi %>%
  dplyr::filter(cell_id %in% raw_cell_ids) %>%
  dplyr::mutate(cell_id=factor(
    as.character(cell_id),
    levels=raw_cell_ids
  )) ->
  df

mat_supernat <- sparseMatrix(
  i=as.integer(df$cell_id),
  j=as.integer(df$Feature),
  x=with(df, both + supernatant_only),
  dims=c(length(levels(df$cell_id)), length(levels(df$Feature))),
  dimnames=list(levels(df$cell_id), levels(df$Feature))
)

mat_pellet <- sparseMatrix(
  i=as.integer(df$cell_id),
  j=as.integer(df$Feature),
  x=with(df, both + pellet_only),
  dims=c(length(levels(df$cell_id)), length(levels(df$Feature))),
  dimnames=list(levels(df$cell_id), levels(df$Feature))
)

mat_both <- sparseMatrix(
  i=as.integer(df$cell_id),
  j=as.integer(df$Feature),
  x=with(df, both + pellet_only + supernatant_only),
  dims=c(length(levels(df$cell_id)), length(levels(df$Feature))),
  dimnames=list(levels(df$cell_id), levels(df$Feature))
)

raw_meta <- data.frame(
  cell_id=raw_cell_ids,
  valid_cell=raw_cell_ids %in% rownames(mdata[['GEX']]$obs)
)

raw_meta$SAMID = stringr::str_match(raw_meta$cell_id, "(SAM.*)_[ACTG]+")[,2]
raw_meta$pass_qc <- pass_qc[as.character(raw_meta$cell_id)]
raw_meta$n_htos <- n_htos_detected[as.character(raw_meta$cell_id)]

stopifnot(raw_meta$cell_id == rownames(mat_both))

sgrna_se <- SummarizedExperiment(
  assays=list(union=mat_both, supernatant=mat_supernat, pellet=mat_pellet),
  rowData=raw_meta
)

stopifnot(rownames(sgrna_se) == rownames(mat_both))
stopifnot(colnames(sgrna_se) == colnames(mat_both))

saveRDS(
  sgrna_se,
  out_se,
  compress=FALSE
)
