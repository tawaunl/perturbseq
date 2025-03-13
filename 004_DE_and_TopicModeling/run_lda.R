library(DataSetDB)
mae <- getDataset("DS000018385", version=3)
sce <- mae[['cellbender']]
rownames(sce) <- scuttle::uniquifyFeatureNames(rownames(sce), rowData(sce)$symbol)

if(!require(fastTopics)){
  install.packages("fastTopics")
  library(fastTopics)
}
print('starting var modeling')
var <- scran::modelGeneVar(sce)
print('finished var modeling')

#chosen <- scran::getTopHVGs(var, n = 1000)
#count.mat <- SummarizedExperiment::assay(sce, "counts")[chosen,]
count.mat <- SummarizedExperiment::assay(sce, "counts")
print('starting LDA model fit')

fit <- fit_topic_model(X = t(as(count.mat, "dgCMatrix")),k = 15)
print('done')

saveRDS(fit, '/gstore/scratch/u/lucast3/ngs5388_collaborative/data/fast_topic_lda_model_fit_all_genes.rds')
