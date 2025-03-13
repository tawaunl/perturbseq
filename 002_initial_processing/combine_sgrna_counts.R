library(magrittr)

#lib_parent_dir <- "/gstore/data/ctgbioinfo/tgi_data_concierge/ngs5388/output"
#libs_csv <- "data/Pilot 9 Trem2 FE perturb - Sample_mapping_for_MThitslims.csv"
#out_path_all <- "scratch/231007_combine_sgrna_counts.csv"
args <- commandArgs(trailingOnly=TRUE)
lib_parent_dir <- args[1]
libs_csv <- args[2]
out_path_all <- args[3]

libs_df <- read.csv(libs_csv)

libs_df %>%
  dplyr::filter(Library.Subtype == 'gRNA V3.1') %>%
  dplyr::rename(portion=Comments) ->
  sgrna_libs_df

sgrna_libs <- with(sgrna_libs_df, paste0(Library.Name, "_", SAMID))
stopifnot(!duplicated(sgrna_libs))
rownames(sgrna_libs_df) <- sgrna_libs

lapply(
  sgrna_libs,
  function(x) {
    readr::read_csv(file.path(
      lib_parent_dir, x, paste0(x, ".stat.csv.gz")
    )) %>%
      dplyr::mutate(lib_sam=x)
  }
) %>%
  do.call(what=rbind) %>%
  dplyr::group_by(lib_sam, Barcode, UMI) %>%
  dplyr::mutate(Total=sum(Count)) %>%
  dplyr::ungroup() %>%
  dplyr::mutate(portion=sgrna_libs_df[lib_sam, "portion"]) %>%
  dplyr::mutate(SAMID=sgrna_libs_df[lib_sam, "SAMID"]) %>%
  dplyr::select(-lib_sam) ->
  sgrna_counts_all

# TODO merge cumulus filtered counts (*.filt.stat.csv.gz)
# TODO merge report files (*.report.txt)

# NOTE: It is recommended to read in the table with stringsAsFactors=TRUE
write.csv(
  sgrna_counts_all,
  out_path_all,
  row.names=FALSE,
  quote=FALSE
)
