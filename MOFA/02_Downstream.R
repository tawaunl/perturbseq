library(ggplot2)
library(MOFA2)
library(dsassembly)
library(SingleCellExperiment)
library(dplyr)
library(tidyverse)
library(Matrix)

reticulate::use_condaenv("~/.conda/envs/mofa_env_gpu/bin/python",required = T)
setwd("/gstore/scratch/u/lucast3/ngs5388_collaborative/MOFA")

# Load Data -----------------
outfile <- paste0(getwd(),"/model_GuideLevelCountswCorrectADT.hdf5")
model <- readRDS(outfile)

mae <- getDataset("DS000018385")
sce <- mae[['cellbender']]


sce_crispr <- altExp(sce, 'CRISPR Guide Capture')
grna_indicator <- as(assay(sce_crispr, 'indicator'),
                     'CsparseMatrix')
geneKO_indicator <- (
  fac2sparse(rowData(sce_crispr)$guide_set_with_ntc_holdouts) %*%
    grna_indicator
) > 0

cite <- altExp(sce, 'Antibody Capture')
# assign KO to metadata
KO <- vector()
for(guide in 1:dim(geneKO_indicator)[1]){
  KO[which(geneKO_indicator[guide,]==1)] <- rownames(geneKO_indicator)[guide]
}

#rename NTCs
KO[which(startsWith(KO,"NTC"))] <- "NTC"
# Add metadata ---------------
sample_metadata <- data.frame(
  sample = samples_names(model)[[1]],
  KO = KO,
  n_genes_prot = cite$n_genes_by_counts,
  total_counts_prot = cite$total_counts
)
sample_metadata <- cbind(sample_metadata,data.frame(colData(sce)))
samples_metadata(model) <- sample_metadata

# Explore Factors and Variance ---------------
cairo_pdf("VarianceDecompositionbyFactor.pdf")
plot_variance_explained(model, x="view", y="factor")+
  theme(axis.text = element_text(size=20),
        axis.text.x = element_text(angle = 45,hjust=.9),
      axis.title = element_text(size=22))
dev.off()

cairo_pdf("CorrelationBetweenFactors.pdf")
plot_factor_cor(model)
dev.off()

cairo_pdf("CountCovariates.pdf")
correlate_factors_with_covariates(model,
                                  covariates = c("n_genes_by_counts","total_counts",
                                                 "n_genes_prot","total_counts_prot"),
                                  fontsize = 20)
dev.off()

cairo_pdf("MetadataCovariates.pdf")
correlate_factors_with_covariates(model,
                                  covariates = c("gate","phase","leiden1"),fontsize = 20)
dev.off()

cairo_pdf("TotalVariance.pdf")
plot_variance_explained(model, plot_total = T)[[2]] +
  theme(axis.text = element_text(size=20),
        axis.title = element_text(size=22))
dev.off()

cairo_pdf("DataOverviewGuideLevel.pdf")
plot_data_overview(model)+
theme(axis.text = element_text(size=20),
      axis.title = element_text(size=22))
dev.off()

# Gate Related Factors -----
plot_factors(model,
             factors = c(3,5,6),
             color_by = "gate"
) + theme(strip.text = element_text(size=18))

# Based on correalition should ignore factor 4
# Factor 1 is interesting related to cell cycle phase and leiden clustering
# Factor 1 Characterization -------------------
plot_factor(model,
            factors = 1,
            color_by = "Factor1"
)

plot_weights(model,
             view = "RNA",
             factor = 1,
             nfeatures = 10,     # Top number of features to highlight
             scale = T           # Scale weights from -1 to 1
)

plot_top_weights(model,
                 view = "RNA",
                 factor = 1,
                 nfeatures = 10,     # Top number of features to highlight
                 scale = T         # Scale weights from -1 to 1
)

plot_factor(model,
            factors = 1,
            color_by = "phase",
            add_violin = TRUE,
            dodge = TRUE
)

# Factor 2 Characterization --------------
# Protein and RNA contributing a good amunt

cairo_pdf("Factor2_TopRNAweights.pdf")
plot_top_weights(model,
                 view = "RNA",
                 factor = 2,
                 nfeatures = 10,     # Top number of features to highlight
                 scale = T         # Scale weights from -1 to 1
)
dev.off()

cairo_pdf("Factor2_TopProteinweights.pdf")
plot_top_weights(model,
                 view = "Protein",
                 factor = 2,
                 nfeatures = 10,     # Top number of features to highlight
                 scale = T         # Scale weights from -1 to 1
)
dev.off()
cairo_pdf("Factor2_TopPerturbweights.pdf")
plot_top_weights(model,
                 view = "Perturbations",
                 factor = 2,
                 nfeatures = 10,     # Top number of features to highlight
                 scale = T         # Scale weights from -1 to 1
)
dev.off()

plot_data_heatmap(model,
                  view = "RNA",
                  factor = 2,
                  features = 25,
                  cluster_rows = FALSE, cluster_cols = FALSE,
                  show_rownames = TRUE, show_colnames = FALSE,
                  scale = "row",denoise = F
)


#Factor 3Characterization ---------------
cairo_pdf("Factor3_TopRNAweights.pdf")
plot_top_weights(model,
                 view = "RNA",
                 factor = 3,
                 nfeatures = 10,     # Top number of features to highlight
                 scale = T         # Scale weights from -1 to 1
)+
  theme(axis.text = element_text(size=20),
        axis.title = element_text(size=22))
dev.off()

cairo_pdf("Factor3_TopProteinweights.pdf")
plot_top_weights(model,
                 view = "Protein",
                 factor = 3,
                 nfeatures = 10,     # Top number of features to highlight
                 scale = T         # Scale weights from -1 to 1
)+
  theme(axis.text = element_text(size=20),
        axis.title = element_text(size=22))
dev.off()

cairo_pdf("Factor3_TopPerturbweights.pdf")
plot_top_weights(model,
                 view = "Perturbations",
                 factor = 3,
                 nfeatures = 10,     # Top number of features to highlight
                 scale = T,sign = "positive"         # Scale weights from -1 to 1
)+
  theme(axis.text = element_text(size=20),
        axis.title = element_text(size=22))
dev.off()

#Factor 5 Characterization ---------------
cairo_pdf("Factor5_TopRNAweights.pdf")
plot_top_weights(model,
                 view = "RNA",
                 factor = 7,
                 nfeatures = 10,     # Top number of features to highlight
                 scale = T         # Scale weights from -1 to 1
)+
  theme(axis.text = element_text(size=20),
        axis.title = element_text(size=22))
dev.off()

cairo_pdf("Factor5_TopProteinweights.pdf")
plot_top_weights(model,
                 view = "Protein",
                 factor = 7,
                 nfeatures = 10,     # Top number of features to highlight
                 scale = T         # Scale weights from -1 to 1
)+
  theme(axis.text = element_text(size=20),
        axis.title = element_text(size=22))
dev.off()

cairo_pdf("Factor5_TopPerturbweights.pdf")
plot_top_weights(model,
                 view = "Perturbations",
                 factor = 7,
                 nfeatures = 10,     # Top number of features to highlight
                 scale = T         # Scale weights from -1 to 1
)+
  theme(axis.text = element_text(size=20),
        axis.title = element_text(size=22))
dev.off()

#Factor 6 Characterization ---------------
## Topweights ----------
cairo_pdf("Factor6_TopRNAweights.pdf")
plot_top_weights(model,
                 view = "RNA",
                 factor = 6,
                 nfeatures = 10,     # Top number of features to highlight
                 scale = T         # Scale weights from -1 to 1
)+
  theme(axis.text = element_text(size=20),
        axis.title = element_text(size=22))
dev.off()

cairo_pdf("Factor6_TopProteinweights.pdf")
plot_top_weights(model,
                 view = "Protein",
                 factor = 6,
                 nfeatures = 10,     # Top number of features to highlight
                 scale = T         # Scale weights from -1 to 1
)+
  theme(axis.text = element_text(size=20),
        axis.title = element_text(size=22))
dev.off()

cairo_pdf("Factor6_TopPerturbweights.pdf")
plot_top_weights(model,
                 view = "Perturbations",
                 factor = 6,
                 nfeatures = 10,     # Top number of features to highlight
                 scale = T,sign = "positive"        # Scale weights from -1 to 1
)+
  theme(axis.text = element_text(size=20),
        axis.title = element_text(size=22))
dev.off()
 ## all weightss ----------
plot_weights(model,
             view="Perturbations",
             factors = 6)
plot_weights_scatter(model,factors = c(3,6),view = "Protein",dot_size=2)

D <- object@dimensions[["D"]][view]
W <- get_weights(object, views=view, factors=factors, as.data.frame = FALSE)
W <- as.data.frame(W); colnames(W) <- c("x","y")
W$view <- view
W$feature <- features_names(object)[[view]]
if (length(color_by)==1 & is.character(color_by)) color_name <- color_by
if (length(shape_by)==1 & is.character(shape_by)) shape_name <- shape_by
color_by <- rep("1",sum(object@dimensions[["D"]][view]))
shape_by <- rep("1",sum(object@dimensions[["D"]][view]))
if (!is(shape_by,"data.frame")) {
  df = data.frame(
    feature = dplyr::filter(model@features_metadata, model@features_metadata$view== "Protein")[,1],
    # sample = unlist(samples_names(object)),
    shape_by = as.factor(shape_by),
    view=view
  )
}
shape_by <- df
if (!is(color_by,"data.frame")) {
  df = data.frame(
    feature = dplyr::filter(model@features_metadata, model@features_metadata$view== "Protein")[,1],
    # sample = unlist(samples_names(object)),
    color_by = as.factor(color_by),
    view=view
  )
}
color_by <- df
# shape_byMerge factor values with group/color/shape information
df <- merge(W, color_by, by=c("feature","view"))
df <- merge(df, shape_by, by=c("feature","view"))

df$shape_by[is.na(df$shape_by)] <- "NA"
df$shape_by <- as.factor(df$shape_by)
if(length(unique(df$color_by)) < 5) df$color_by <- as.factor(df$color_by)
df$x <- df$x/max(abs(df$x))
df$y <- df$y/max(abs(df$y))

p <- ggplot(df, aes(x=.data$x, y=.data$y,labels=feature)) +
  geom_point(aes(color = .data$color_by, shape = .data$shape_by), size=3) +
  labs(x=factors[1], y=factors[2]) +
  geom_segment(x=min(df$x,na.rm=TRUE), xend=max(df$x,na.rm=TRUE), y=0, yend=0, linewidth=0.25, color="orange") +
  geom_segment(y=min(df$y,na.rm=TRUE), yend=max(df$y,na.rm=TRUE), x=0, xend=0, linewidth=0.25, color="orange") +
  theme_classic() +
  theme(
    axis.text = element_text(size=rel(1), color="black"),
    axis.title = element_text(size=rel(1.3), color="black"),
    axis.ticks = element_line(color="black")
  )+
  coord_cartesian(xlim=c(-1,1), ylim=c(-1,1))+
  guides(color="none") + scale_color_manual(values="black")+
  theme(legend.position = "none") +
  ggrepel::geom_text_repel(label=df$feature,size=7,max.overlaps = 4)

p


plot_weights_scatter(model,factors = c(3,6),view = "RNA",dot_size=2)
plot_weights_scatter(model,factors = c(3,6),view = "Perturbations",dot_size=2)
object <- model
factors = c(3,6)
view = "Perturbations"
shape_by <- NULL
color_by <- NULL
D <- object@dimensions[["D"]][view]
W <- get_weights(object, views=view, factors=factors, as.data.frame = FALSE)
W <- as.data.frame(W); colnames(W) <- c("x","y")
W$view <- view
W$feature <- features_names(object)[[view]]
if (length(color_by)==1 & is.character(color_by)) color_name <- color_by
if (length(shape_by)==1 & is.character(shape_by)) shape_name <- shape_by
color_by <- rep("1",sum(object@dimensions[["D"]][view]))
shape_by <- rep("1",sum(object@dimensions[["D"]][view]))
if (!is(shape_by,"data.frame")) {
  df = data.frame(
    feature = dplyr::filter(model@features_metadata, model@features_metadata$view== "Perturbations")[,1],
    # sample = unlist(samples_names(object)),
    shape_by = as.factor(shape_by),
    view=view
  )
}
shape_by <- df
if (!is(color_by,"data.frame")) {
  df = data.frame(
    feature = dplyr::filter(model@features_metadata, model@features_metadata$view== "Perturbations")[,1],
    # sample = unlist(samples_names(object)),
    color_by = as.factor(color_by),
    view=view
  )
}
color_by <- df
# shape_byMerge factor values with group/color/shape information
df <- merge(W, color_by, by=c("feature","view"))
df <- merge(df, shape_by, by=c("feature","view"))

df$shape_by[is.na(df$shape_by)] <- "NA"
df$shape_by <- as.factor(df$shape_by)
if(length(unique(df$color_by)) < 5) df$color_by <- as.factor(df$color_by)
df$x <- df$x/max(abs(df$x))
df$y <- df$y/max(abs(df$y))
labels <- df[which(abs(df$x)>0.25),]
labels <- rbind(labels,df[which(abs(df$y)>0.25),])
p <- ggplot(df, aes(x=.data$x, y=.data$y,labels=feature)) +
  ggrastr::geom_point_rast(aes(color = .data$color_by, shape = .data$shape_by), size=3) +
  labs(x=factors[1], y=factors[2]) +
  geom_segment(x=min(df$x,na.rm=TRUE), xend=max(df$x,na.rm=TRUE), y=0, yend=0, linewidth=0.25, color="orange") +
  geom_segment(y=min(df$y,na.rm=TRUE), yend=max(df$y,na.rm=TRUE), x=0, xend=0, linewidth=0.25, color="orange") +
  theme_classic() +
  theme(
    axis.text = element_text(size=rel(1), color="black"),
    axis.title = element_text(size=rel(1.3), color="black"),
    axis.ticks = element_line(color="black")
  )+ggrepel::geom_text_repel(data = labels,label=labels$feature,size=7,max.overlaps = 10) +
  coord_cartesian(xlim=c(-1,1), ylim=c(-1,1))+
  guides(color="none") + scale_color_manual(values="black")+
  theme(legend.position = "none")

p

labels <- df[which(df$x>0.25),]
labels <- rbind(labels,df[which(df$y>0.25),])

factors <- 1:get_dimensions(model)[["K"]]
factors <- factors[!factors%in%c(4,7)]
