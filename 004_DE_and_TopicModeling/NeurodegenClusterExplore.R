# Neurodegen cluster explore

library(DataSetDB)
mae <- getDataset("DS000016134")
sce <- mae[['gene']]

lda_res <- readRDS('/gstore/scratch/u/lucast3/TREM2perturb/ngs5388_collaborative/fast_topic_lda_model_fit_all_genes.rds')

# Cluster marker analysis ------------------
library(scran.chan)
counts <- counts(sce)
rownames(counts) <- make.unique(rowData(sce)$symbol)
counts <- initializeSparseMatrix(counts)
lognorm <- logNormCounts.chan(counts,batch = sce$assignment)
markers <- scoreMarkers.chan(lognorm,
                             groups=sce$leiden1, batch=sce$assignment,lfc=0)

library(edgeR)
# Get top 5 markers from each cluster
topmarkers <- c()
for (cluster in 1:length(markers$statistics)) {
  clustermarkers <- data.frame(markers$statistics[[cluster]])
  clustermarkers <- clustermarkers[order(clustermarkers$logFC,decreasing = TRUE),]
  top <- rownames(clustermarkers)[1:5]
  topmarkers <- c(topmarkers,top)

}


library(RColorBrewer)
library(viridis)
abundances <- table(sce$leiden1,sce$SAMID)
abundances <- unclass(abundances)

extra.info <- colData(sce)[match(colnames(abundances), sce$SAMID),]
d <- DGEList(abundances, samples=extra.info)
d = calcNormFactors(d)
d= estimateCommonDisp(d, verbose=TRUE)



norm_counts <- as.data.frame(t(d$counts))
colnames(norm_counts) <- paste("Cluster ", colnames(norm_counts), sep="")
norm_counts <- norm_counts %>% mutate(SAMID = rownames(.))
coldata <- as.data.frame(colData(sce))
coldata_short <- coldata %>% dplyr::select(SAMID,gate) %>% unique


df_long_final <- left_join(norm_counts,coldata_short, by="SAMID")

percentages <- data.frame(matrix(nrow = dim(df_long_final)[1],ncol = dim(df_long_final)[2]))
sums <- data.frame(clustersums=rowSums(df_long_final[,1:dim(norm_counts)[2]-1]),ident=df_long_final$SAMID)

for (clust in 1:length(levels(factor(sce$leiden1)))) {
  for (sample in 1:length(df_long_final[,clust])) {
    percent <- df_long_final[sample,clust]/sums[sample,1]*100
    percentages[sample,clust]<- percent
  }
}
cols2replace <- dim(norm_counts)[2]:(dim(norm_counts)[2]+dim(coldata_short)[2]-1)
percentages[,cols2replace]<- df_long_final[,cols2replace]


colnames(percentages) <-  colnames(df_long_final)
df_long <- percentages %>%
  pivot_longer(c(c(paste0("Cluster ",0:(dim(norm_counts)[2]-2)))))

colnames(df_long) <- c("Sample","Pool","Cluster","Percent")

level_order <- c("Trem2Hi","Trem2Lo","All")
df_long$Pool <- factor(df_long$Pool,levels = level_order)
df_long$Cluster <- factor(df_long$Cluster,levels =c(paste0("Cluster ",0:(dim(norm_counts)[2]-2))))



ggplot(df_long, aes(x=Pool, y=Percent,fill=Pool,colour = Pool)) + theme_classic() +
  geom_jitter(color="darkgrey",width = 0.4,size=3,alpha=0.6)  +
  geom_boxplot(outlier.shape=NA,alpha=0.5) +
  facet_wrap(~factor(Cluster),scales = "free", ncol=4) + xlab("Pool") +
  theme(axis.text.y = element_text(size=20)) +
  ylab("Cellularity Percent") +
  theme(strip.text.x = element_text(size = 20,face = "bold"),
        axis.title=element_text(size=16,face="bold"),
        legend.text = element_text(size=14),
        legend.title = element_text(size=16),
        axis.title.x = element_blank(),
        axis.text.x = element_blank())
# DE on cluster 6 vs all -------------------------------------
library(scran)
library(edgeR)
summed <- aggregateAcrossCells(sce,
                               id=colData(sce)[,c("SAMID","leiden1","gate")])

summed.filt <- summed[,summed$ncells >= 10]
current <- summed.filt
y <- DGEList(counts(current), samples=colData(current))
keep <- filterByExpr(y, group=current$leiden1)
y <- y[keep,]
y <- calcNormFactors(y)

design <- model.matrix(~0+leiden1+gate, data=y$samples)
colnames(design) <- paste0("Cluster_",levels(y$samples$leiden1))
my.contrasts <- makeContrasts((Cluster_0+Cluster_1 +Cluster_2+Cluster_3+
                                 Cluster_4 + Cluster_5+ Cluster_7 + Cluster_8+
                                 Cluster_9 + Cluster_10+ Cluster_11)/11-Cluster_6,
                              levels=colnames(design))
y <- estimateDisp(y, design)
fit <- glmQLFit(y, design, robust=TRUE)



lrt <- glmQLFTest(fit, contrast=my.contrasts)
#res <- glmQLFTest(fit, coef=ncol(design),contrast = my.contrasts)

res <- topTags(lrt, n = Inf)$table

res <- res %>% mutate(sig = ifelse((PValue<=0.05 | PValue <= 0.05) & abs(logFC) >= 1, "yes", "no"))
x <- match(rownames(res),rownames(rowData(sce)))
syms <- rowData(sce)$gene_name[x]
id <- rownames(rowData(sce))[x]
res$symbol <- syms
res$ID <- id

source("~/scHelpers.R")
rownames(res) <- res$symbol
res$logFC <- res$logFC*-1
edgeRVolcano(res,plot.title = "Cluster6 vs ALL",labels = rownames(res)[which(res$sig=="yes")])
saveRDS(res,"/gstore/scratch/u/lucast3/TREM2perturb/ngs5388_collaborative/Tawaun_analysis/Cluster6vsAll_DE.rds")

topFeatures <- rownames(res)[which(res$sig=="yes" & res$logFC > 1.5)]

library(scater)
rownames(sce) <- make.unique(rowData(sce)$gene_name)
plotDots(summed,features = topFeatures[1:50],group="leiden1",exprs_values="counts",
         center=T,scale=T)
saveRDS(summed,"/gstore/scratch/u/lucast3/TREM2perturb/ngs5388_collaborative/Tawaun_analysis/PuesdobulkClusters.rds")
## Pathway analysis on Top genes -------------------------------
library(clusterProfiler)
library(enrichplot)
res <- readRDS("/gstore/scratch/u/lucast3/TREM2perturb/ngs5388_collaborative/Tawaun_analysis/Cluster6vsAll_DE.rds")


res <- res[order(res$logFC,decreasing = TRUE),]
diffGenes <- res$logFC
names(diffGenes) <- res$ID

# convert gene ID to ENTREZID
gene.df <- bitr( names(diffGenes), fromType = "ENSEMBL",
                 toType = c("SYMBOL", "ENTREZID"),
                 OrgDb = org.Mm.eg.db::org.Mm.eg.db)
x <- match(gene.df$ENSEMBL,names(diffGenes))
gseaVec <- diffGenes[x]
names(gseaVec) <- gene.df$SYMBOL
gseaVec <- gseaVec[order(gseaVec,decreasing = T)]
# GSEA analysis, using fgsea under hood
ego3 <- gseGO(geneList     = gseaVec,
              OrgDb        = org.Mm.eg.db::org.Mm.eg.db,
              ont          = "BP",
              keyType = "SYMBOL",
              minGSSize    = 100,
              maxGSSize    = 500,
              eps          = 0,
              pvalueCutoff = 0.05,
              verbose      = FALSE)
paths <- dotplot(filter(ego3,NES >0),x="NES",color="pvalue",decreasing=TRUE)+
  ggtitle("Cluster6 GO:Enriched Pathways") +
  theme(axis.title.x = element_text(size = 14,face="bold"),
        axis.text.x = element_text(angle = 90)) +
  theme(plot.title = element_text(hjust = 0.5,size=20,face=20))

ego3 <- pairwise_termsim(filter(ego3,NES>0))
treeplot(filter(ego3,NES>0))

# Explore Topics ---------------------------------
library(grid)
library(gridExtra)
library(ggpubr)
do_gsea <- function(topic, i_, Pathways, minSize = 15, maxSize = 500){
  dat <- topic[,i_]
  #print(dat)
  fgseaRes <- fgsea(pathways = Pathways,
                    stats    = dat,
                    minSize  = minSize,
                    maxSize  = maxSize)
  return(fgseaRes)
}
# join loadings to sce object
topic_range <- colnames(lda_res$L)
for (topic in topic_range) {
  sce[[paste0("lda_",topic)]] <- lda_res$L[,topic]
}

# Prepare list for scores
n_genes = 25
n_features = 200
topic_score_dict = list()
plot_score_dict = list()

for (topic in topic_range) {
  topic_scores = lda_res$F[,topic]
  topic_scores <- sort(topic_scores,decreasing = T)
  topic_dict = list()
  for (i in 1:n_features) {
    gene <- names(topic_scores)[i]
    score = as.numeric(topic_scores[i])
    topic_dict[gene] = score
  }
  topic_score_dict[[topic]] = topic_dict

  plot_series = lapply(topic_dict,sort,decreasing=TRUE)
  plot_score_dict[[topic]] = plot_series

}

topic_plots <- list()
x=1
for(topic in c("k7","k8")){
  df <- data.frame(gene=names(topic_score_dict[[topic]]),
                   score=unlist(topic_score_dict[[topic]]))
  topic_plots[[x]] <- ggplot(df[1:20,],aes(x=reorder(gene,score),y=score)) +
    geom_bar(stat="identity") +
    coord_flip() +
    scale_y_continuous(expand = c(0, 0),
                       labels = function(x) format(x, scientific = TRUE),
                       breaks = range(df$score)) +
    xlab("Gene") + ggtitle(paste0("Topic",x)) +
    theme_light()+ theme(plot.title = element_text(hjust = 0.5))
  x=x+1
}

do.call("grid.arrange",list(grobs=topic_plots,top=textGrob("Gene scores for each topic")))

topic_plots_umap <- list()
x=1
for(topic in topic_range){

  topic_plots_umap[[x]] <- plotReducedDim(sce,
                                          colour_by = paste0("lda_",topic),
                                          dimred = "X_wnn_umap") +
    ggtitle(paste0("Topic ",x)) +
    theme_void() +
    theme(plot.title = element_text(hjust = 0.5,face = "bold"), legend.position = "none")
  x=x+1
}
### Topic 14 Focus -----------------------------
df <- data.frame(gene=names(topic_score_dict[["k5"]]),
                 score=unlist(topic_score_dict[["k5"]]))
ggplot(df[1:25,],aes(x=reorder(gene,score),y=score)) +
  geom_bar(stat="identity") +
  coord_flip() +
  scale_y_continuous(expand = c(0, 0),
                     labels = function(x) format(x, scientific = TRUE),
                     breaks = range(df$score)) +
  xlab("Gene") + ggtitle(paste0("Topic",5)) +
  theme_light()+ theme(plot.title = element_text(hjust = 0.5))

gsea_res <- lapply(1:15, function(x){do_gsea(fit$F, x, go_list)})

gsea_14 <- gsea_res[[9]]
topPathwaysUp <- gsea_14[ES > 0][head(order(pval), n=15), pathway]

plotGseaTable(go_list[topPathwaysUp], fit$F[,9], gsea_14,
              gseaParam=0.5)

collapsedPathways <- collapsePathways(gsea_14[order(pval)][padj < 0.01],
                                      go_list, lda_res$F[,14])
mainPathways <- gsea_14[pathway %in% collapsedPathways$mainPathways][
  order(-NES), pathway]
plotGseaTable(go_list[mainPathways], lda_res$F[,14], gsea_14,
              gseaParam = 0.5)
#subcluster ----------------
library(Seurat)
source("~/scHelpers.R")
cluster6 <- sce[,sce$leiden1==6]
saveRDS(cluster6,"/gstore/scratch/u/lucast3/TREM2perturb/ngs5388_collaborative/Tawaun_analysis/Cluster6Subset_sce.rds")
counts(cluster6) <- as(counts(cluster6),"dgCMatrix")
data <- CreateSeuratObject(counts = counts(cluster6),meta.data = data.frame(colData(cluster6)))
data <- runSeurat(data)
scvi_coords <- get_scvi_coords(data,data$seurat_clusters)
ggplot(data=scvi_coords , aes(x=UMAP1, y=UMAP2, colour  = seurat_clusters)) +
  geom_point(size=2,alpha = 0.6) +theme_classic()
ggplot(data=scvi_coords , aes(x=UMAP1, y=UMAP2, colour  = gate)) +
  geom_point(size=2,alpha = 0.6) +theme_classic()

## Diff. Abundance --------------
abundances <- table(data$seurat_clusters,data$SAMID)
abundances <- unclass(abundances)

extra.info <- data@meta.data[match(colnames(abundances), data$SAMID),]
d <- DGEList(abundances, samples=extra.info)
d = calcNormFactors(d)
d= estimateCommonDisp(d, verbose=TRUE)



norm_counts <- as.data.frame(t(d$counts))
colnames(norm_counts) <- paste("Cluster ", colnames(norm_counts), sep="")
norm_counts <- norm_counts %>% mutate(SAMID = rownames(.))
coldata <- as.data.frame(data@meta.data)
coldata_short <- coldata %>% dplyr::select(SAMID,gate) %>% unique


df_long_final <- left_join(norm_counts,coldata_short, by="SAMID")

percentages <- data.frame(matrix(nrow = dim(df_long_final)[1],ncol = dim(df_long_final)[2]))
sums <- data.frame(clustersums=rowSums(df_long_final[,1:dim(norm_counts)[2]-1]),ident=df_long_final$SAMID)

for (clust in 1:length(levels(factor(data$seurat_clusters)))) {
  for (sample in 1:length(df_long_final[,clust])) {
    percent <- df_long_final[sample,clust]/sums[sample,1]*100
    percentages[sample,clust]<- percent
  }
}
cols2replace <- dim(norm_counts)[2]:(dim(norm_counts)[2]+dim(coldata_short)[2]-1)
percentages[,cols2replace]<- df_long_final[,cols2replace]


colnames(percentages) <-  colnames(df_long_final)
df_long <- percentages %>%
  pivot_longer(c(c(paste0("Cluster ",0:(dim(norm_counts)[2]-2)))))

colnames(df_long) <- c("Sample","Pool","Cluster","Percent")

level_order <- c("Trem2Hi","Trem2Lo","All")
df_long$Pool <- factor(df_long$Pool,levels = level_order)
df_long$Cluster <- factor(df_long$Cluster,levels =c(paste0("Cluster ",0:(dim(norm_counts)[2]-2))))

ggplot(df_long, aes(x=Pool, y=Percent,fill=Pool,colour = Pool)) + theme_classic() +
  geom_jitter(color="darkgrey",width = 0.4,size=3,alpha=0.6)  +
  geom_boxplot(outlier.shape=NA,alpha=0.5) +
  facet_wrap(~factor(Cluster),scales = "free", ncol=3) + xlab("Pool") +
  theme(axis.text.y = element_text(size=20)) +
  ylab("Cellularity Percent") +
  theme(strip.text.x = element_text(size = 20,face = "bold"),
        axis.title=element_text(size=16,face="bold"),
        legend.text = element_text(size=14),
        legend.title = element_text(size=16),
        axis.title.x = element_blank(),
        axis.text.x = element_blank())
## Topic Modeling --------------------
rownames(cluster6) <- scuttle::uniquifyFeatureNames(rownames(cluster6), rowData(cluster6)$symbol)

print('starting var modeling')
var <- scran::modelGeneVar(cluster6)
print('finished var modeling')

chosen <- scran::getTopHVGs(var, n = 1000)
count.mat <- SummarizedExperiment::assay(cluster6, "counts")[chosen,]
#count.mat <- SummarizedExperiment::assay(sce, "counts")
print('starting LDA model fit')

fit <- fit_topic_model(X = t(as(count.mat, "dgCMatrix")),k = 15)
print('done')

saveRDS(fit, '/gstore/scratch/u/lucast3/TREM2perturb/ngs5388_collaborative/Cluster6_fast_topic_lda_model_fit.rds')
### Explore topics ---------------
# join loadings to sce object
topic_range <- colnames(fit$L)
for (topic in topic_range) {
  cluster6[[paste0("lda_",topic)]] <- fit$L[,topic]
  data[[paste0("lda_",topic)]] <- fit$L[,topic]
}

# Prepare list for scores
n_genes = 25
n_features = 200
topic_score_dict = list()
plot_score_dict = list()

for (topic in topic_range) {
  topic_scores = fit$F[,topic]
  topic_scores <- sort(topic_scores,decreasing = T)
  topic_dict = list()
  for (i in 1:n_features) {
    gene <- names(topic_scores)[i]
    score = as.numeric(topic_scores[i])
    topic_dict[gene] = score
  }
  topic_score_dict[[topic]] = topic_dict

  plot_series = lapply(topic_dict,sort,decreasing=TRUE)
  plot_score_dict[[topic]] = plot_series

}

topic_plots <- list()
x=1
for(topic in topic_range){
  df <- data.frame(gene=names(topic_score_dict[[topic]]),
                   score=unlist(topic_score_dict[[topic]]))
  topic_plots[[x]] <- ggplot(df[1:20,],aes(x=reorder(gene,score),y=score)) +
    geom_bar(stat="identity") +
    coord_flip() +
    scale_y_continuous(expand = c(0, 0),
                       labels = function(x) format(x, scientific = TRUE),
                       breaks = range(df$score)) +
    xlab("Gene") + ggtitle(paste0("Topic",x)) +
    theme_light()+ theme(plot.title = element_text(hjust = 0.5))
  x=x+1
}

do.call("grid.arrange",list(grobs=topic_plots,top=textGrob("Gene scores for each topic")))

topic_plots_umap <- list()
x=1
for(topic in topic_range){

  topic_plots_umap[[x]] <- FeaturePlot(data,
                                       paste0("lda_",topic),
                                       pt.size = 1.5,order=T) +
    ggtitle(paste0("Topic ",x)) +
    theme_void() +
    theme(plot.title = element_text(hjust = 0.5,face = "bold"),
          legend.position = "none")
  x=x+1
}

do.call("grid.arrange",
        list(grobs=topic_plots_umap,
             top=textGrob(expression(bold(underline("UMAP embedding for each topic"))))))

## look at kO
plt_kos <- c('Urod', 'Furin',
             'Eef2', 'Lamtor1', 'Lamtor2', 'Lamtor3', 'Rraga',
             'Cdkn1a', 'Fgfr1op', 'Ppp2r3c', 'Sec62', 'Sec63',
             'Trem2', 'NTC')

sce_crispr <- altExp(cluster6, 'CRISPR Guide Capture')

# matrix of guide-level calls
grna_calls <- assay(sce_crispr, 'indicator')

# combine guides to gene level
modmat <- fac2sparse(rowData(sce_crispr)$guide_set)
geneKO_calls <- (modmat %*% grna_calls) > 0
geneKO_calls <- geneKO_calls[-c(1:4),]

ko <- cbind(
  data@meta.data[, c("seurat_clusters","gate"), drop=FALSE],
  data@reductions$umap@cell.embeddings,
  t(as.matrix(geneKO_calls))
) %>%
  as.data.frame() %>%
  tidyr::pivot_longer(dplyr::all_of(rownames(geneKO_calls)),
                      names_to='KO', values_to='indicator') %>%
  dplyr::arrange(indicator) %>%
  dplyr::mutate(KO=factor(KO, levels=rownames(geneKO_calls))) %>%
  dplyr::mutate(gate=dplyr::if_else(
    indicator, gate, NA_character_
  ))
test <- ko %>%
  group_by(seurat_clusters, KO,indicator) %>%
  summarise(Count = n())
colnames(test)[1] <- "Cluster"
test1 <- test[test$indicator=="TRUE",]

test1 <- test1[,-3]
test1$Count <- as.numeric(test$Count)
library(dplyr)
library(tidyverse)
dat <-data.frame(pivot_wider(test1,names_from = KO,values_from = Count))
dat$Cluster <- paste0("Cluster_",dat$Cluster)
rownames(dat) <-dat$Cluster
dat<- as.matrix(dat[,-1])
dat[which(is.na(dat))] <-as.numeric(0)
dat1 <- dat[,-455] # dont plot NTC
library(pheatmap)
pheatmap(dat1,scale = "column",)

scaled <- scale(dat[,-14])

# DiffExpr. Trem2 hi vs low --------------
# DE on cluster 6 vs all -------------------------------------
library(scran)
library(edgeR)
summed <- aggregateAcrossCells(sce,
                               id=colData(sce)[,c("SAMID","gate")])

summed.filt <- summed[,summed$ncells >= 10]
current <- summed.filt
y <- DGEList(counts(current), samples=colData(current))
keep <- filterByExpr(y, group=current$gate)
y <- y[keep,]
y <- calcNormFactors(y)

design <- model.matrix(~0+gate, data=y$samples)
colnames(design) <- paste0("Cluster_",levels(y$samples$gate))
my.contrasts <- makeContrasts(Cluster_Trem2Hi-Cluster_Trem2Lo,
                              levels=colnames(design))
y <- estimateDisp(y, design)
fit <- glmQLFit(y, design, robust=TRUE)

lrt <- glmQLFTest(fit, contrast=my.contrasts)
#res <- glmQLFTest(fit, coef=ncol(design),contrast = my.contrasts)

res <- topTags(lrt, n = Inf)$table

res <- res %>% mutate(sig = ifelse((PValue<=0.05 | PValue <= 0.05) & abs(logFC) >= 1, "yes", "no"))
x <- match(rownames(res),rownames(rowData(sce)))
syms <- rowData(sce)$gene_name[x]
id <- rownames(rowData(sce))[x]
res$symbol <- syms
res$ID <- id

source("~/scHelpers.R")
rownames(res) <- make.unique(res$symbol)
cluster6_res <- readRDS("/gstore/scratch/u/lucast3/TREM2perturb/ngs5388_collaborative/Tawaun_analysis/Cluster6vsAll_DE.rds")

library(gg4way)

data <- list(Trem2HivsLow=res,
             Cluster6vsAll= cluster6_res)

gg4way(DGEdata = data,
       x = "Trem2HivsLow",
       y = "Cluster6vsAll", sep = "vs",FDR="PValue",
       logFC = "logFC",logFCcutoff = 1, label=T, textSize =14,
       colorVector = c("grey80", "darkgreen", "mediumblue", "red")) +
  xlab(expression(atop(
    paste("Higher in Trem2Low" %<->% "Higher in TREM2Hi"),
    paste("Trem2 Hi Macrophges vs Trem2 Low  LogFC")))) +
  ylab(expression(atop(
    paste("Cluster 6 Macrophages vs All other Clusters LogFC"),
    paste("Higher in Other Clusters" %<->% "Higher in Cluster 6")))) +
  theme(axis.title = element_text(size=22),
        axis.text = element_text(size=20))


