library(magrittr)

library(Matrix)

library(SummarizedExperiment)
library(SingleCellExperiment)
library(MultiAssayExperiment)

library(dsassembly)

print("Read the SCE")

initial_processing_config <- yaml::read_yaml('../002_initial_processing/config.yaml')
sce <- readRDS(file.path(initial_processing_config$scratch_dir, 'convert_mdata_to_sce.Rds'))

print("Load old version of data from DSDB")

mae_old <- getDataset("DS000016134", version=1)
sce_old <- mae_old[['gene']]

print("Fixup rowdata of new SCE")

stopifnot(nrow(sce_old) == nrow(sce))
stopifnot(!duplicated(rownames(sce)))
stopifnot(setequal(rowData(sce_old)$old_name, rownames(sce)))

sce <- sce[rowData(sce_old)$old_name,]
rownames(sce) <- rownames(sce_old)

for (k in c("gene_id", "symbol", "gene_name")) {
  rowData(sce)[,k] <- rowData(sce_old)[,k]
}

print("Annotate RNA")

sce %<>% dsassembly::annotateExperiment(
  title = "Filtered barcode-gene counts",
  description = "Barcode-gene counts from cumulus with ambient background subtraction by cellbender.",
  organism = "Mus musculus",
  annotation = metadata(sce_old)$.internal$DataSetDB$summarized_experiment$annotation$id,
  technology = list(name="scRNA-seq", details="10X Genomics"),
  sources = list(list(name="FireDB", id="FRS19421"))
)

print("Annotate HTO")

altExp(sce, 'Multiplexing Capture') %<>% dsassembly::annotateExperiment(
  title = "Barcode-HTO counts",
  description = "Barcode-HTO counts from cumulus.",
  organism = "Mus musculus",
  technology = list(name="CellHashing", details="10X Genomics"),
  annotation=NULL,
  namespace=list(),
  sources = list(list(name="FireDB", id="FRS19423"))
)

print("Annotate ADT")

altExp(sce, 'Antibody Capture') %<>% dsassembly::annotateExperiment(
  title = "Barcode-ADT counts",
  description = "Barcode-ADT counts from cumulus.",
  organism = "Mus musculus",
  technology = list(name="CITE-seq", details="10X Genomics"),
  annotation=NULL,
  namespace=list(),
  sources = list(list(name="FireDB", id="FRS19422"))
)

print("Annotate gRNA")

altExp(sce, 'CRISPR Guide Capture') %<>% dsassembly::annotateExperiment(
  title = "Barcode-guide assignments",
  description = "Barcode-guide counts from cumulus.",
  organism = "Mus musculus",
  technology = list(name="gRNA-seq", details="10X Genomics"),
  annotation='GMTY236:Trem_pilot7_sgRNA_lib@REVISION-1',
  sources = list(list(name="FireDB", id="FRS19424"), list(name="FireDB", id="FRS19425"))
)

print("Make sample mapping")

mapping <- data.frame(
  assay = "cellbender",
  primary = colData(sce)$SAMID,
  colname = colnames(sce)
)

print("Make MAE colData")

colData(sce) %>%
  as.data.frame() %>%
  dplyr::select(SAMID, sample_name) %>%
  dplyr::distinct() %>%
  `rownames<-`(NULL) %>%
  dplyr::arrange(SAMID) %>%
  dplyr::mutate(SPECIES='Mus musculus') ->
  mae_coldat

rownames(mae_coldat) <- mae_coldat$SAMID

print("Create MultiAssayExperiment")

mae <- MultiAssayExperiment(
  experiments=list(cellbender=sce),
  colData=mae_coldat,
  sampleMap=mapping
)

print("Annotate MAE")

mae <- dsassembly::annotateDataset(
  mae,
  title = "NGS5388: Flow enriched perturbseq for Trem2 in mouse BMDMs (cellbender-cleaned version)",
  description = "Trem2 hypofunctional mutations increase risk of developing Alzheimers Disease 3fold. We are using guideRNA perturbations for genes relevant to Trem2 biology and AD to identify novel regulators of Trem2 surface localization and function. Briefly, we collect bone marrow from mouse femurs and tibias and differentiate them into macrophages. We transduce a 5000 guide lentiviral library and nucleofect the cells with Cas9. Cells are hashed (6HTOs) and stained with Biolegend Total SeqA CITE-Seq panel as well as Trem2 CITE-Seq compatible antibody to further divide cells into high and low Trem2 categories. We collect cells representative of all the population, cells expressing High levels of Trem2 and cells expressing low levels of Trem2. This version of the dataset uses cellbender ambient-subtracted gene expression counts; an older version of this dataset, without cellbender, is DS000016134.",
  authors = "kammj2"
)

print("Upload to datasetdb")

# in case any h5 matrices are higher than the 1Gb default limit
options(ArtifactDB.upload.size.limit = 5)

dsassembly::saveDataset(mae)

# Succesfully uploaded as:
# DS000018385
