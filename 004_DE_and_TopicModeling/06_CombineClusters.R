library(DataSetDB)

mae <- getDataset("DS000018385")
sce <- mae[['cellbender']]


library(scran)
library(BiocParallel)

#out_leiden_1 <- findMarkers(sce, groups=sce$leiden1)
#out_leiden_2 <- findMarkers(sce, groups=sce$leiden2,BPPARAM=bp)


#save(out_leiden_1, out_leiden_2, file = '/gstore/scratch/u/lucast3/ngs5388_collaborative/scratch/single_cell_markers.rdata')
load('/gstore/scratch/u/lucast3/ngs5388_collaborative/scratch/single_cell_markers.rdata')
rownames(sce) <- scuttle::uniquifyFeatureNames(rownames(sce), rowData(sce)$symbol)

out_leiden_1_changed <- lapply(out_leiden_1, function(x){
  x$gene_id <- rownames(x)

  x <- merge(x, rowData(sce), by = 'gene_id')
  rownames(x) <- x$symbol
  x <- x[order(x$Top),]
  return(x)
})

is_distinct <- function(cluster_name, marker_res, .se, cluster = 'leiden1', top_n = 30, assay_to_use = "logcounts"){
  markers <- rownames(head(subset(marker_res, marker_res$summary.logFC>0), n= top_n))
  count.mat <- assay(.se, 'logcounts')[markers,]
  average.mat <- data.frame(gene = markers)
  for(u_cluster in unique(colData(.se)[[cluster]])){
    cells <- rownames(colData(.se)[colData(.se)[[cluster]] == u_cluster,])
    col_to_add <- as.data.frame(rowMeans(count.mat[,cells]))
    names(col_to_add) <- c(u_cluster)
    average.mat<-cbind(average.mat, col_to_add)
  }
  average.mat$gene <- NULL
  effects <- average.mat - rowMeans(average.mat, na.rm=TRUE)
  pheatmap(effects)
  dist <- dist(t(effects) , diag=TRUE)
  hc <- hclust(dist)
  cut_avg <- cutree(hc, k = 2)
  unique_cluster_ <- names(cut_avg[cut_avg == names(table(cut_avg)[table(cut_avg) == 1])])
  if(length(unique_cluster_)== 0){
    print("there are no unique cluster")
    return(FALSE)
  }
  if(cluster_name == unique_cluster_){
    return(TRUE)
  } else {
    print(paste0("this is unique?? ", unique_cluster_))
    return(FALSE)
  }
}
is_distinct("0", out_leiden_1_changed[['0']], sce, top_n = 10)
is_distinct("1", out_leiden_1_changed[['1']], sce, top_n = 10) #TRUE
is_distinct("2", out_leiden_1_changed[['2']], sce, top_n = 10)
is_distinct("3", out_leiden_1_changed[['3']], sce, top_n = 10)
is_distinct("4", out_leiden_1_changed[['4']], sce, top_n = 10)
is_distinct("5", out_leiden_1_changed[['5']], sce, top_n = 10)
is_distinct("6", out_leiden_1_changed[['6']], sce, top_n = 10)
is_distinct("7", out_leiden_1_changed[['7']], sce, top_n = 10)
is_distinct("8", out_leiden_1_changed[['8']], sce, top_n = 10) #TRUE
is_distinct("9", out_leiden_1_changed[['9']], sce, top_n = 10)
is_distinct("10", out_leiden_1_changed[['10']], sce, top_n = 10)
is_distinct("11", out_leiden_1_changed[['11']], sce, top_n = 10)
is_distinct("12", out_leiden_1_changed[['12']], sce, top_n = 10)



# this thing is taking a while so lets try and do something else.
# look at correlation of FC??
library(pheatmap)
dat <- data.frame(place_holder =out_leiden_1_changed[[1]]$gene_id)
for(cluster in 0:13){
  data <- out_leiden_1_changed[[as.character(cluster)]]
  data <- data[order(data$gene_id),]
  data_ <- data.frame("a" = data$summary.logFC)
  names(data_) <- as.character(cluster)
  dat <- cbind(dat, data_)
}
dat$place_holder <- NULL
cor_dat <- cor(as.matrix(dat))
#cor_dat[cor_dat<0.7] = 0
png("heatmapClusters.R",height = 2000,width=2000, res=150)
pheatmap(cor_dat,display_numbers = T)
dev.off()
saveRDS(cor_dat,"/gstore/scratch/u/lucast3/ngs5388_collaborative/corDat.rds")

# OK I think we can merge 7,1,2??
# maybe there is a good cutoff for correaltion ? 0.6

mae <- getDataset("DS000018385")
sce <- mae[['cellbender']]
sce$round1_merge <- sce$leiden1
sce$round1_merge <- gsub("2", "1", sce$round1_merge)
sce$round1_merge <- gsub("7", "1", sce$round1_merge)
#finally combine 8 & 9
sce$round1_merge <- gsub("9", "8", sce$round1_merge)

# Rerun Correlation ?
outround1 <- readRDS("/gstore/scratch/u/lucast3/ngs5388_collaborative/round1.rds")

round1_change <- lapply(outround1, function(x){
  x$gene_id <- rownames(x)
  x <- x[order(x$Top),]
  return(x)
})
dat <- data.frame(place_holder =round1_change[[1]]$gene_id)


for(cluster in names(round1_change)){
  data <- round1_change[[as.character(cluster)]]
  data <- data[order(data$gene_id),]
  data_ <- data.frame("a" = data$summary.logFC)
  names(data_) <- as.character(cluster)
  dat <- cbind(dat, data_)
}
dat$place_holder <- NULL
cor_datr1 <- cor(as.matrix(dat))
pheatmap(cor_datr1,display_numbers = T)

# not too convinced if we need to combine any other clusters

saveRDS(sce,"./scratch/WorkingSCE.rds")
