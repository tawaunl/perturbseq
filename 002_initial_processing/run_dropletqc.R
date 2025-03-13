library(magrittr)

# TODO: Just use velocyto instead. DropletQC is not packaged well (has
# to be manually installed), and their docs even recommend velocyto
# (https://powellgenomicslab.github.io/DropletQC/articles/DropletQC.html)
library(DropletQC)

#lib_csv <- "data/Pilot 9 Trem2 FE perturb - Sample_mapping_for_MThitslims.csv"
#parent_dir <- "/gstore/data/ctgbioinfo/tgi_data_concierge/ngs5388/output"
#out_rds <- "scratch/231005_dropletqc.Rds"
args <- commandArgs(trailingOnly=TRUE)
lib_csv <- args[1]
parent_dir <- args[2]
out_rds <- args[3]

ncores <- as.integer(Sys.getenv()[["SLURM_CPUS_ON_NODE"]])

lib_df <- read.csv(lib_csv)
lib_df %<>% dplyr::filter(Library.Subtype == 'Expression V3.1')

stopifnot(!duplicated(lib_df$SAMID))
stopifnot(!duplicated(lib_df$Sample.Name))

lib_sam <- with(lib_df, paste0(Library.Name, '_', SAMID))

gex_folders <- file.path(parent_dir, lib_sam)
names(gex_folders) <- lib_df$SAMID

# Use verbose=FALSE to prevent dropletqc from spamming error log with
# its progress bar

res <- lapply(
  gex_folders,
  function(f) nuclear_fraction_tags(
    f, cores=ncores,
    verbose=FALSE
  )
)

saveRDS(res, out_rds, compress = FALSE)
