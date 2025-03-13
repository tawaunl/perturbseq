library(magrittr)
library(Matrix)

library(reticulate)
use_python(python=Sys.which("python"))

muon <- import("muon")

args <- commandArgs(trailingOnly=TRUE)
in_mdata <- args[1]
out_mdata <- args[2]
out_confusion_csv <- args[3]
out_best_params_csv <- args[4]
denom_pseudo_max <- as.integer(args[5])

#mdata <- muon$read("scratch/refine_dbl.h5mu")
mdata <- muon$read(in_mdata)

mdata[['GEX']]$obs %>%
  dplyr::mutate(
    total_grna_umi=mdata[['sgRNA']]$obs$n_umi,
    first_grna_umi=mdata[['sgRNA']]$obs$max_umi_cnt,
    second_grna_umi=mdata[['sgRNA']]$obs$second_umi_cnt
  ) %>%
  dplyr::filter(
    qc_cluster_dbl_refined %in% c('qc_clust1', 'hto_doublet_11')
  ) %>%
  .[, c('sample_name', 'demux_type', 'total_grna_umi',
        'first_grna_umi', 'second_grna_umi')] %>%
  dplyr::mutate(
    is_dbl=demux_type == 'doublet'
  ) ->
  df_mod_dbl_sgrna

get_confusion <- function(n_cutoff, frac_cutoff, denom_pseudo=0) {
  df_mod_dbl_sgrna %>%
    dplyr::mutate(
      pred_demux=factor(dplyr::case_when(
        first_grna_umi < n_cutoff | first_grna_umi / (total_grna_umi+denom_pseudo) < frac_cutoff ~ 'unassigned',
        second_grna_umi < n_cutoff | second_grna_umi / (total_grna_umi+denom_pseudo) < frac_cutoff ~ 'singlet',
        TRUE ~ 'doublet'
      ), levels=c('singlet', 'doublet', 'unassigned'))
    ) %>%
    .[, c('demux_type', 'pred_demux')] %>%
    table()
}

expand.grid(
  n_cutoff=1:5,
  #frac_cutoff=seq(0, .5, .05),
  frac_cutoff=seq(0, 50, 5) / 100,
  denom_pseudo=0:denom_pseudo_max
) %>%
  as.data.frame() ->
  params_grid

apply(params_grid, 1, function(x) {
  get_confusion(
    x['n_cutoff'], x['frac_cutoff'], x['denom_pseudo']
  ) %>%
    as.data.frame() %>%
    dplyr::mutate(
      n_cutoff=x['n_cutoff'],
      frac_cutoff=x['frac_cutoff'],
      denom_pseudo=x['denom_pseudo']
    )
}) %>%
  do.call(what=rbind) ->
  df_confusion

df_confusion %>%
  tidyr::pivot_wider(
    names_from=c(demux_type, pred_demux),
    values_from=Freq
  ) %>%
  dplyr::mutate(dbl_frac_correct = doublet_doublet / (doublet_doublet +
                                                        doublet_singlet +
                                                        doublet_unassigned)) %>%
  dplyr::filter(dbl_frac_correct >= .5) %>%
  dplyr::filter(n_cutoff >= 2) %>%
  dplyr::arrange(desc(singlet_singlet)) %>%
  .[1,] ->
  best_params

grna_cnt <- mdata[['sgRNA']]$X

grna_calls <- (grna_cnt >= best_params[['n_cutoff']]) & (
  (Diagonal(x=1/(rowSums(grna_cnt) + best_params[['denom_pseudo']]))
    %*% grna_cnt) >= best_params[['frac_cutoff']]
)

grna_calls %<>% drop0()

## convert to numeric (do we need to do this?)
#grna_calls <- grna_calls + 0.0

py_run_string("r.mdata['sgRNA'].raw = r.mdata['sgRNA']")
py_run_string("r.mdata['sgRNA'].X = r.grna_calls")

mdata$write(out_mdata)

write.csv(
  df_confusion,
  out_confusion_csv,
  row.names=F
)

write.csv(
  best_params,
  out_best_params_csv,
  row.names=F
)
