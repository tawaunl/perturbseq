library(magrittr)

library(data.table)
library(dtplyr)
library(dplyr)

source("/gstore/scratch/u/lucast3/ngs5388_collaborative/030_rna_diff_expr/paths.sh")

list.files(file.path(ANALYSIS_SCRATCH, "020_run_glmgampoi_singlecell_genelvl"),
           full.names=T) %>%
  lapply(fread) %>%
  rbindlist() ->
  res_glmgp

out_path <- file.path(ANALYSIS_RESULTS, "021_merge_genelvl_de_results.csv")
fwrite(res_glmgp, out_path)
