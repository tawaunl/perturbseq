library(magrittr)
library(Matrix)

library(reticulate)
use_python(python=Sys.which("python"))

muon <- import("muon")

args <- commandArgs(trailingOnly=TRUE)
in_mdata <- args[1]
out_prefix <- args[2]
out_cutoffs_csv <- args[3]

#mdata <- muon$read("scratch/add_sgrna_mudata.h5mu", backed=TRUE)
mdata <- muon$read(in_mdata)

X_pca <- mdata[['GEX']]$obsm[['X_pca']]

is_dbl <- mdata[['GEX']]$obs$demux_type == 'doublet'

res_lda <- MASS::lda(
  X_pca[!is_dbl,],
  mdata[['GEX']]$obs$qc_cluster_raw[!is_dbl]
)

ld_pred <- predict(res_lda, X_pca)$x[,1]

qc_clust_raw <- mdata[['GEX']]$obs$qc_cluster_raw
# not sure if lda actually enforces the direction, but we assume it
# below, so we assert it. If direction turns out to be random, and
# group1 can be higher than group2, then multiply ld_pred by -1 in
# that case.
stopifnot(
  median(ld_pred[qc_clust_raw == 1]) <
    median(ld_pred[qc_clust_raw == 2])
)

# cutoffs based on per-sample fractions

unique(mdata[['GEX']]$obs$sample_name) %>%
  as.character() ->
  sample_names

sample_names %>%
  lapply(function(x) {
    keep <- mdata[['GEX']]$obs$sample_name == x

    n <- table(qc_clust_raw[!is_dbl & keep])

    p11 <- n[1]^2 / sum(n)^2
    p22 <- n[2]^2 / sum(n)^2

    data.frame(
      lower=quantile(ld_pred[is_dbl & keep], p11),
      upper=quantile(ld_pred[is_dbl & keep], 1-p22),
      lower_quantile=p11,
      upper_quantile=1-p22
    )
  }) %>%
  do.call(what=rbind) ->
  sample_cutoffs

rownames(sample_cutoffs) <- sample_names

mdata[['GEX']]$obs[, 'qc_cluster_dbl_refined'] <- with(
  mdata[['GEX']]$obs,
  dplyr::case_when(
    qc_cluster != 'hto_doublet' ~ qc_cluster,
    ld_pred < sample_cutoffs[as.character(sample_name), 'lower'] ~ 'hto_doublet_11',
    ld_pred > sample_cutoffs[as.character(sample_name), 'upper'] ~ 'hto_doublet_22',
    TRUE ~ 'hto_doublet_12'
  )
)

mdata[['GEX']]$obs[, 'qc_dbl_lda_pred'] <- ld_pred

# save the updated mudata
mdata$write(paste0(out_prefix, '.h5mu'))

saveRDS(
  res_lda,
  paste0(out_prefix, '_lda.Rds'),
  compress = FALSE
)

sample_cutoffs %>%
  as.data.frame() %>%
  tibble::rownames_to_column("sample_name") %>%
  write.csv(
    out_cutoffs_csv,
    row.names=F
  )
