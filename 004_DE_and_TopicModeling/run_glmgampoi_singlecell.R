# Originally copied from:
# https://code.roche.com/kammj2/ngs5388_trem2feps_pilot9/-/blob/e8df89648828908c24282188b33c51ac5455a33f/src/run_glmgampoi_singlecell.R

library(magrittr)
library(Matrix)

library(SummarizedExperiment)
library(SingleCellExperiment)
library(MultiAssayExperiment)

library(glmGamPoi)

library(BiocParallel)

library(optparse)

arguments <- parse_args(
  OptionParser(
    usage = "%prog [options] sce_rds out_csv",
    option_list = list(
      make_option("--n_subproc", type="integer", default=1),
      make_option("--array_jobs", type="integer", default=1),
      make_option("--array_idx", type="integer", default=1),
      make_option("--ridge_prior_sigma", type="double", default=0),
      make_option("--gene_min_cnt", type="integer", default=100),
      make_option("--use_shared_overdispersion",
                  action="store_true", default=FALSE,
                  help="Estimate overdispersion on the control cells only."),
      make_option("--no_shrink_overdisp",
                  action="store_true", default=FALSE,
                  help="Don't apply shrinkage to overdispersion estimates"),
      make_option("--test",
                  action="store_true", default=FALSE,
                  help="Run a smaller test job"),
      make_option("--guide_lvl",
                  action="store_true", default=FALSE,
                  help="Test at guide level (instead of gene level).")
    )
  ),
  positional_arguments=2
)

args <- arguments$args
opt <- arguments$options

sce_rds <- args[1]
out_csv <- args[2]

use_shared_overdispersion <- opt$use_shared_overdispersion
no_shrink_overdisp <- opt$no_shrink_overdisp
test_guide_lvl <- opt$guide_lvl
test_subset <- opt$test
n_subproc <- opt$n_subproc
array_idx <- opt$array_idx
array_jobs <- opt$array_jobs
ridge_prior_sigma <- opt$ridge_prior_sigma
gene_min_cnt <- opt$gene_min_cnt

sce <- readRDS(sce_rds)

rna_coldat <- colData(sce)

X_cnt_raw <- counts(sce)

keep_genes <- rowSums(X_cnt_raw) >= gene_min_cnt
print("Number of genes to keep")
summary(keep_genes)

X_cnt_raw %<>% .[keep_genes,]

grna_indicator <- assay(altExp(sce, 'CRISPR Guide Capture'), 'indicator')

stopifnot(colnames(grna_indicator) == colnames(X_cnt_raw))

crispr_rowdat <- rowData(altExp(sce, 'CRISPR Guide Capture'))

ctrl_guides <- crispr_rowdat$guide_set_with_ntc_holdouts == 'NTC_1_100'
ctrl_cells <- colSums(grna_indicator[ctrl_guides,]) >= 1

print("Number of control cells")
summary(ctrl_cells)

stopifnot(rownames(rna_coldat) == colnames(X_cnt_raw))
sample_design_mat <- model.matrix(~0+sample_name, rna_coldat)
cells_per_guide <- rowSums(grna_indicator)

if (use_shared_overdispersion) {
  res_ctrl_glmgampoi <- glmGamPoi::glm_gp(
    X_cnt_raw[, ctrl_cells], design = sample_design_mat[ctrl_cells,],
    on_disk=FALSE
    #on_disk=TRUE
  )
  disp_est <- res_ctrl_glmgampoi$overdispersions
}

ko_mat <- grna_indicator

if (!test_guide_lvl) {
  modmat <- fac2sparse(crispr_rowdat$guide_set_with_ntc_holdouts)
  ko_mat <- (modmat %*% ko_mat) > 0
  ko_mat %<>% .[rownames(.) != 'NTC_1_100',]

  if (test_subset) {
    ko_mat <- ko_mat[c('Trem2', 'NTC_101_104', 'Ran', '1700020N01Rik'),]
  }
} else {
  ko_mat %<>% .[!(rownames(.) %in% paste0('NTC_', 1:100)),]

  if (test_subset) {
    ko_mat <- ko_mat[c(paste0('Trem2_',1:4),
                       paste0('NTC_', 101:104),
                       'NTC_189',
                       'Ran_2',
                       paste0('1700020N01Rik_', 1:4)),]
  }
}

# remove empty KOs
ko_mat %<>% .[rowSums(.) > 0,]

if (array_jobs > 1) {
  keep_ko <- ((1:nrow(ko_mat)) - 1) %% array_jobs == (array_idx - 1)
  ko_mat <- ko_mat[keep_ko,]
}

run_glmgp <- function(curr_ko) {
  #curr_ko <- 'Ran'
  #curr_ko <- 'Trem2'
  #curr_ko <- 'NTC_101_104'

  ko_cells <- ko_mat[curr_ko,] >= 1
  union_cells <- ctrl_cells | ko_cells

  design_mat <- cbind(
    sample_design_mat[union_cells,],
    as.integer(ko_mat[curr_ko, union_cells])
  )
  colnames(design_mat) <- make.names(c(colnames(sample_design_mat), curr_ko))

  curr_cnt <- X_cnt_raw[, union_cells]

  # get the ridge penalty
  if (ridge_prior_sigma != 0) {
    ridge_pen <- c(
      rep(0, ncol(sample_design_mat)),
      # equivalent to Gaussian prior on log2FC whose SD is
      # ridge_prior_sigma
      log(2) / ridge_prior_sigma / sqrt(nrow(design_mat))
    )
  } else {
    ridge_pen <- 0
  }

  if (use_shared_overdispersion) {
    res_glmgp <- glm_gp(
      curr_cnt, design=design_mat,
      overdispersion = disp_est,
      overdispersion_shrinkage = !no_shrink_overdisp,
      ridge_penalty = ridge_pen,
      on_disk=FALSE
    )
  } else {
    res_glmgp <- glm_gp(
      curr_cnt, design=design_mat,
      overdispersion_shrinkage = !no_shrink_overdisp,
      ridge_penalty = ridge_pen,
      on_disk=FALSE
    )
  }

  test_de(res_glmgp,
          contrast=make.names(curr_ko),
          max_lfc=Inf) %>%
    as.data.frame() %>%
    dplyr::mutate(KO=curr_ko) ->
    ret

  # get standard error
  idmat <- diag(nrow=ncol(res_glmgp$Beta))
  colnames(idmat) <- colnames(res_glmgp$Beta)
  rownames(idmat) <- colnames(res_glmgp$Beta)

  pred <- predict(
    res_glmgp,
    se.fit=TRUE,
    newdata = idmat
  )

  stopifnot(ret$name == rownames(pred$se.fit))

  ret$se <- pred$se.fit[,make.names(curr_ko)] / log(2)

  ret
}

sprintf("Running DE on %d KOs", nrow(ko_mat))

if (n_subproc > 1) {
  # Limit num cpus per process with RhpcBLASctl
  # https://stackoverflow.com/questions/45794290/in-r-how-to-control-multi-threading-in-blas-parallel-matrix-product
  library(RhpcBLASctl)
  blas_set_num_threads(1)
  omp_set_num_threads(1)

  system.time(
    combined_results <- bplapply(
      rownames(ko_mat), run_glmgp,
      BPPARAM=MulticoreParam(n_subproc)
    )
  )
} else {
  system.time(
    combined_results <- lapply(
      1:nrow(ko_mat),
      function(i) {
        x <- rownames(ko_mat)[i]

        # TODO Add verbosity flag to control this
        message(sprintf(
          "[%s] %d / %d (%s)",
          format(Sys.time(), "%Y-%m-%d %H:%M:%S"),
          i, nrow(ko_mat), x
        ))

        run_glmgp(x)
      }
    )
  )
}


combined_results <- do.call(rbind, combined_results)

# mkdir -p
dir.create(dirname(out_csv), recursive=TRUE)

write.csv(
  combined_results,
  out_csv,
  row.names=F
)
