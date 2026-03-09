library(SingleCellExperiment)
library(Signac)
library(dplyr)
library(ggplot2)
library(scDblFinder)
library(scuttle)
library(scran)
library(edgeR)
library(ggrepel)
library(EDASeq)
library(muscat)

multiome_sce <- readRDS("~/multiome_sce.rds")
multiome_sce
#multiome_sham <- multiome_sce[,multiome_sce$condition == "SHAM"]

#PCA per vedere se ci sono effetti di batch
pb_mult <- aggregateAcrossCells(multiome_sce, use.assay.type = "counts", 
                                id= multiome_sce$dataset)
pb_mult <- logNormCounts(pb_mult)
pb_mult <- runPCA(pb_mult, BSPARAM = BiocSingular::ExactParam(), ncomponents=2)
plotPCA(pb_mult, color_by = "dataset", point_size=8)


#Provo a rifare i DEG con RUV
library(RUVSeq)
library(scSEG)
get_seg_controls <- function() {
  data(segList)
  segList$mouse$mouse_scSEG
}
controls <- get_seg_controls()

pb_mult <- aggregateAcrossCells(multiome_sham, use.assay.type = "counts", 
                                id=DataFrame(label = multiome_sham$subclusters,
                                             sample = multiome_sham$dataset))

pb_mult$group <- as.factor(pb_mult$subclusters)
y <- DGEList(counts(pb_mult), samples = as.data.frame(colData(pb_mult)))
keep <- filterByExpr(y, group = y$samples$group)
y <- y[keep,, keep.lib.sizes=FALSE]
y <- calcNormFactors(y)
design <- model.matrix(~ 0 + group, data = y$samples)
y <- estimateDisp(y, design)
fit <- glmFit(y, design)

# RUVg factor
r <- RUVg(y$counts, control_genes, k = k)
y$samples$W <- r$W[,1, drop=TRUE]
design <- model.matrix(~ 0 + group + W, data = y$samples) 
colnames(design) <- levels(y$samples$group) #?
y <- estimateDisp(y, design)
fit <- glmFit(y, design)



