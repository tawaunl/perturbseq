library(fastTopics)
library(DataSetDB)
mae <- getDataset("DS000018385")
sce <- mae[['cellbender']]

lda_res <- readRDS('/gstore/scratch/u/lucast3/ngs5388_collaborative/data/fast_topic_lda_model_fit_all_genes.rds')


library(dplyr)

# GO lists
library(biomaRt)
ensembl <- useEnsembl(biomart = "genes", dataset = "mmusculus_gene_ensembl")
data <- getBM(attributes = c('external_gene_name', 'name_1006'),
              filters = 'external_gene_name',
              values = rownames(lda_res$F),
              mart = ensembl)
data <- data[data$name_1006 != "",]

go_list <- list()
for(go_term in unique(data$name_1006)){
  d <- data[data$name_1006 == go_term,]
  go_list <- append(go_list, list(d$external_gene_name))
}
names(go_list)<- unique(data$name_1006)
read.csv('/gstore/scratch/u/lucast3/ngs5388_collaborative/data/friedman_modules.csv') %>%
  apply(2, function(x) stringr::str_to_title(x[x != ""])) ->
  gene_modules

go_list <- append(go_list, gene_modules)

for_enricher <- data.frame()
for(term in names(go_list)){
  for_enricher <- rbind(for_enricher, data.frame(term = term , genes = go_list[[term]]))
}
library(clusterProfiler)
library(fgsea)
do_gsea <- function(topic, i_, Pathways, minSize = 15, maxSize = 500){
  dat <- topic[,i_]
  print(dat)
  fgseaRes <- fgsea(pathways = Pathways,
                    stats    = dat,
                    minSize  = minSize,
                    maxSize  = maxSize)
  return(fgseaRes)
}
do_enrich_go <- function(topics, i_, go_enrich, top_n_genes = 200){
  df <- topics[,i_]
  df <- df[order(-df)]
  results <- enricher(names(head(df, top_n_genes)),
                      universe = names(df),
                      gson = NULL,
                      TERM2GENE= go_enrich)
  return(results@result)
}
enrich_res <- lapply(1:15, function(x){do_enrich_go(lda_res$F, x, for_enricher, top_n_genes = 50)})
gsea_res <- lapply(1:15, function(x){do_gsea(lda_res$F, x, go_list)})

library(igraph)
source("/gstore/scratch/u/lucast3/ngs5388_collaborative/scratch/makeTopics.R")
gsea_res_ <- lapply(gsea_res, function(x){return(as.data.frame(x))})
network <- make_topics_network(go_list_ = enrich_res,
                               col_name_terms = 'ID',
                               col_name_pvalue = 'p.adjust',
                               pvalue_cutoff = 0.05, overlap_lim = 3)
V(network)$size  <- V(network)$size *5
V(network)$name[V(network)$size  < 5] <- NA
V(network)$size[V(network)$size  < 5] <- 0.0000001
V(network)$color <- unlist(lapply(V(network)$type, function(x){if(x=="topics"){return('salmon')}else{return("light blue")}}))
l_IP <- layout_with_fr(network)
plot(network, layout=l_IP, vertex.label.family = "Arial", vertex.label.cex=0.8)

# cell cycle
plot(lda_res$F[,2], lda_res$F[,1])
library(ggplot2)
# Dam Cells
a <- data.frame(lda_res$F)
a <- a[order(-a$k1),]
a$gene <- rownames(a)
a$gene <- factor(a$gene, levels=a$gene)
ggplot(a[1:20,], aes(gene, k1)) + geom_bar(stat = 'identity')
View(enrich_res[[11]])

lda_res$L
coldata <- as.data.frame(colData(sce))
sum(rownames(coldata) != rownames(lda_res$L))
coldata <- cbind(coldata, lda_res$L)
library(tidyr)
coldata_ <- coldata %>%
  pivot_longer(
    cols = starts_with("k"),
    names_to = "topic",
    values_to = "score",
    values_drop_na = TRUE
  )
ggplot(coldata_, aes(x = round1_merge, y = score)) + geom_boxplot() +
  facet_wrap(~topic) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
ggplot(coldata, aes(x = round1_merge, y = k15)) + geom_boxplot()
ggplot(coldata, aes(x = round1_merge, y = k11)) + geom_boxplot()


rownames(sce) <- scuttle::uniquifyFeatureNames(rownames(sce), rowData(sce)$symbol)

plot(lda_res$F['Lgals3',][order(-lda_res$F['Lgals3',])])
#plotExpression(sce, x ="round1_merge", features = c("Ctsd"))

coldata$topic_cell_cycle <- coldata$k3 + coldata$k6 #blue
coldata$topic_neurodegen <- coldata$k7 + coldata$k8
coldata$topic_Spp1 <- coldata$k2 #red
coldata$topic_chemokine <- coldata$k5 #magenta
coldata$topic_mito <-  coldata$k15
coldata$topic_ribo <-  coldata$k11
coldata$topic_mac <- coldata$k13 + coldata$k14
coldata$topic_Ctsd <- coldata$k1
coldata$topic_Lgals3_Prdx1 <- coldata$k4
coldata$topic_unknown <- coldata$k9 + coldata$k10 + coldata$k12


coldata$id <- rownames(coldata)
coldata_anno_topic <- coldata %>%
  pivot_longer(
    cols = starts_with('topic'),
    names_to = "topic",
    values_to = "score",
    values_drop_na = TRUE
  )



ggplot(coldata_anno_topic, aes(x = round1_merge, y = score)) + geom_boxplot() +
  facet_wrap(~topic) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

lda_res_to_test <- readRDS('/gstore/scratch/u/seos6/ngs5388_collaborative/fast_topic_lda_model_fit.rds')
lda_res_to_test$L <- as.matrix(coldata[startsWith(names(coldata),'topic')])

f_tb <- as.data.frame(lda_res_to_test$F)

f_tb$topic_cell_cycle <- (f_tb$k3 + f_tb$k6)/2 #blue
f_tb$topic_neurodegen <- (f_tb$k7 + f_tb$k8)/2
f_tb$topic_Spp1 <- f_tb$k2 #red
f_tb$topic_chemokine <- f_tb$k5 #magenta
f_tb$topic_mito <-  f_tb$k15
f_tb$topic_ribo <-  f_tb$k11
f_tb$topic_mac <- (f_tb$k13 + f_tb$k14)/2
f_tb$topic_Ctsd <- f_tb$k1
f_tb$topic_Lgals3_Prdx1 <- f_tb$k4
f_tb$topic_unknown <- (f_tb$k9 + f_tb$k10 + f_tb$k12)/3

lda_res_to_test$F <- as.matrix(f_tb[startsWith(names(f_tb),'topic')])



p1a <- structure_plot(lda_res_to_test)
p1b <- structure_plot(lda_res_to_test, grouping = sce$round1_merge, gap = 25)

#lets update the mae
mae <- getDataset("DS000016134")
sce <- mae[['gene']]
sum(rownames(colData(sce)) != rownames(coldata))
ph <- colData(sce)
ph <- cbind(ph, coldata[startsWith(names(coldata),'topic')])
colData(sce) <- ph
mae[['gene']] <- sce
updateDataset(mae)


write.csv(f_tb[startsWith(names(f_tb),'topic')], file = '/gstore/scratch/u/seos6/bmdm_trem2_topics.csv')

