library(magrittr)

library(data.table)
library(dtplyr)
library(dplyr)

library(Matrix)

library(MultiAssayExperiment)
library(SingleCellExperiment)

library(dsassembly)

# git clone https://code.roche.com/kammj2/xshrink.git
# devtools::install('/path/to/xshrink')
library(xshrink)

library(BiocParallel)

if ("SLURM_CPUS_ON_NODE" %in% names(Sys.getenv())) {
  ncores <- as.integer(Sys.getenv()[["SLURM_CPUS_ON_NODE"]])
  BPPARAM <- MulticoreParam(ncores)
} else {
  BPPARAM <- MulticoreParam(4)
}

source("/gstore/scratch/u/lucast3/ngs5388_collaborative/030_rna_diff_expr/paths.sh")

## Load the SCE object and extract some objects for plotting
mae <- getDataset("DS000018385")
sce <- mae[['cellbender']]
rownames(sce) <- make.unique(rowData(sce)$symbol)

sce_crispr <- altExp(sce, 'CRISPR Guide Capture')

# make sure cells are in the same order across altExps
stopifnot(colnames(sce) == colnames(sce_crispr))

# grna by cell, 1=cell has grna, 0=else
grna_indicator <- as(assay(sce_crispr, 'indicator'),
                     'CsparseMatrix')

# indicator matrix pooling all guides within a gene KO
geneKO_indicator <- (
  fac2sparse(rowData(sce_crispr)$guide_set_with_ntc_holdouts) %*%
    grna_indicator
) > 0

## Load the DE results

res_glmgp <- fread(file.path(ANALYSIS_RESULTS, "021_merge_genelvl_de_results.csv"))

res_glmgp2 <- lazy_dt(res_glmgp)

## Extract matrices of p-values, FDR, LFCs, SEs

res_glmgp2 %>%
  dplyr::select(KO, name, adj_pval) %>%
  as.data.frame() %>%
  tidyr::pivot_wider(names_from=KO, values_from=adj_pval) %>%
  tibble::column_to_rownames('name') %>%
  as.matrix() ->
  mat_qvals

res_glmgp2 %>%
  dplyr::select(KO, name, lfc) %>%
  as.data.frame() %>%
  tidyr::pivot_wider(names_from=KO, values_from=lfc) %>%
  tibble::column_to_rownames('name') %>%
  as.matrix() ->
  mat_lfc

res_glmgp2 %>%
  dplyr::select(KO, name, pval) %>%
  as.data.frame() %>%
  tidyr::pivot_wider(names_from=KO, values_from=pval) %>%
  tibble::column_to_rownames('name') %>%
  as.matrix() ->
  mat_pval

res_glmgp2 %>%
  dplyr::select(KO, name, se) %>%
  as.data.frame() %>%
  tidyr::pivot_wider(names_from=KO, values_from=se) %>%
  tibble::column_to_rownames('name') %>%
  as.matrix() ->
  mat_se

## Shrink the LFCs

res_shrunk <- xshrink(
  mat_lfc, mat_se
)

mat_shrunk <- res_shrunk$shrunken

## Compute ICA

ica_shrunk <- ica::icafast(
  t(mat_shrunk), nc=8,
  center=FALSE
)

## Plot heatmap

# plotting KOs with at least 10 DE genes
plt_kos <- colSums(mat_qvals <= .05) >= 10
plt_kos %<>% {names(.[.])}

# Plot the top 20 genes in each ICA component
plt_genes <- apply(
  lda_res$F, 2,
  function(x) {
    x <- sort(abs(x))
    tail(names(x), n=20)
    #x <- sort(x)
    #c(head(names(x), n=10), tail(names(x), n=10))
  }
)

plt_genes <- unique(unlist(as.data.frame(plt_genes)))

# LFCs to plot on heatmap
mat_plt <- mat_shrunk[plt_genes,][,plt_kos]

# Put a star or X on significant entries of the heatmap
mat_star <- matrix("", nrow=nrow(mat_plt), ncol=ncol(mat_plt),
                   dimnames = dimnames(mat_plt))
mat_star[mat_qvals[plt_genes,][,plt_kos] <= .05] <- 'x'

# Row annotations (ICA loadings)
lda_res$F %>%
  as.data.frame() %>%
  .[rownames(mat_plt),] %>%
  `colnames<-`(sprintf('LDA_%02d', 1:ncol(.))) %>%
  as.data.frame() ->
  anno_row

# Column annotations

# number of cells in each KO and gate
n_per_gate <- geneKO_indicator[plt_kos,] %*% t(fac2sparse(colData(sce)$gate))
# for each KO, the fraction of cells in each gate
frac_per_gate <- sweep(n_per_gate, 1, rowSums(n_per_gate), '/')
frac_per_gate <- as.data.frame(as.matrix(frac_per_gate))

anno_col <- cbind(num_de_log10=log10(colSums(mat_qvals <= .05)[plt_kos]+1),
                  frac_per_gate)

p <- ComplexHeatmap::Heatmap(
  mat_plt,
  clustering_method_rows = 'ward.D2',
  clustering_method_columns = 'ward.D2',
  row_names_gp = grid::gpar(fontsize=8),
  column_names_gp = grid::gpar(fontsize=10),
  col = circlize::colorRamp2(c(-2,0,2),
                             c('blue','white','red')),
  cell_fun = function(j, i, x, y, w, h, col) { # add text to each grid
    grid::grid.text(mat_star[i, j], x, y,
                    gp=grid::gpar(fontsize=8))
  },
  top_annotation = ComplexHeatmap::HeatmapAnnotation(
    df=anno_col,
    col=c(list(num_de_log10=circlize::colorRamp2(seq(0, 3, length.out=11),
                                                 pals::brewer.greys(11))),
          apply(frac_per_gate, 2,
                function(x) circlize::colorRamp2(seq(0, 1, length.out=10),
                                                 pals::parula(10)),
                simplify=FALSE)),
    show_legend=c(TRUE, TRUE, FALSE, FALSE),
    #show_legend=FALSE,
    annotation_legend_param = apply(frac_per_gate, 2,
                                    function(x) list(title='Fraction'),
                                    simplify=FALSE)
  ),
  left_annotation = ComplexHeatmap::HeatmapAnnotation(
    df=anno_row,
    col=apply(anno_row, 2,
              function(x) circlize::colorRamp2(seq(-.15,.15, length.out=11),
                                               pals::coolwarm(11)),
              simplify=FALSE),
    #show_legend=FALSE,
    show_legend=c(TRUE, rep(FALSE, ncol(anno_row)-1)),
    annotation_legend_param = apply(anno_row, 2,
                                    function(x) list(title='Loading'),
                                    simplify=FALSE),
    which='row'
  ),
  heatmap_legend_param = list(title = 'LFC')
)

png("/gstore/scratch/u/lucast3/ngs5388_collaborative/scratch/results/heatmap_shrunk_lfc_with_LDA.png",
    width=1300, height=1200)
ComplexHeatmap::draw(p, merge_legend=TRUE)
dev.off()
