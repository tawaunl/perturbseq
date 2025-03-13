library(magrittr)
library(Matrix)

library(reticulate)
use_python(python=Sys.which("python"))

muon <- import("muon")

#in_mdata = "scratch/assign_grnas.h5mu"
#keep_adt_thresh <- 0.95

args <- commandArgs(trailingOnly=TRUE)
in_mdata <- args[1]
out_mdata <- args[2]
keep_adt_thresh <- as.numeric(args[3])

mdata <- muon$read(in_mdata)

# reticulate has trouble with sparse integer matrices. So convert the
# matrix in Python before moving to R
#X_adt_full_raw <- mdata[['ADT']]$raw$X
print('Converting/moving ADT counts')
py_run_string("X = r.mdata['ADT'].raw.X.todense()")
X_adt_full_raw <- py$X

keep_cells <- with(mdata[['GEX']]$obs, qc_cluster=='qc_clust1')

keep_adt <- colMeans(X_adt_full_raw[keep_cells,] > 0) >= keep_adt_thresh

print("ADTs above cutoff:")
summary(keep_adt)

X_sub <- X_adt_full_raw[keep_cells,][,keep_adt]
X_sub <- log(X_sub + 1)

res_medpolish <- medpolish(
  X_sub, eps=1e-6
)

log(X_adt_full_raw[,keep_adt] + 1) %>%
  sweep(2, res_medpolish$col + res_medpolish$overall, '-') %>%
  apply(1, median) ->
  log_size_factors

stopifnot(all.equal(
  log_size_factors[keep_cells],
  res_medpolish$row +
    # this is 0 iff the last step of polish was row-wise
    apply(res_medpolish$residual, 1, median)
))

X_medpolish <- log1p(sweep(
  X_adt_full_raw, 1, exp(log_size_factors), '/'
))

py_run_string("r.mdata['ADT'].X = r.X_medpolish")

stopifnot(all.equal(X_medpolish, mdata[['ADT']]$X))

# TODO save the size factors

mdata$write(out_mdata)
