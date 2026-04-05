#fermiamoci un attimo a guardare l'espressione dei vari tipi cellulari

#16bin ---------------------------------------------------------------------------------------------------------------
bin16_list_ann <- readRDS("~/bin16_list_ann.rds")

#mettiamo assieme gli sham e partiamo da là
bin16_list_ann <- lapply(bin16_list_ann, function(spe) {
   keep <- !is.na(spe$new_cell_type) & colSums(counts(spe)) > 0
   spe[,keep]
})
for (nm in names(bin16_list_ann)) {
  spe <- bin16_list_ann[[nm]]
  colData(spe)$sample_id <- rep(nm, ncol(spe))
  bin16_list_ann[[nm]] <- spe
}
bin16_list_ann$c26SMAD23_b2$to_discard <- NULL
bin16_list_ann$c26murf1_b2$to_discard <- NULL
combined_spe <- do.call(cbind, bin16_list_ann)
sce_sham <- combined_spe[,combined_spe$tissue_type == "sham"]

#average expression for each gene in each celltype
library(scuttle)
ids <- paste(
  sce_sham$sample_id,
  sce_sham$new_cell_type,
  sep = "_"
)
avg_expr <- aggregateAcrossCells(
  sce_sham,
  ids = ids,
  use.assay.type = "logcounts"
)

#number of detected genes
library(Matrix)
counts_mat <- assay(sce_sham, "counts") > 0
genes_per_cell <- colSums(counts_mat)
sce_sham$detected_genes <- genes_per_cell
tapply(sce_sham$detected_genes,
       sce_sham$new_cell_type,
       mean)
#for example
sce_sham[,sce_sham$cell_type == "MuSC"]$detected_genes

# most expressed genes for each cell type
# top_genes <- apply(assay(avg_expr), 2, function(x) {
#   names(sort(x, decreasing = TRUE))[1:10]
# })
# top_genes

#but maybe if I divided for libsize
mat <- assay(avg_expr)
# scaling per gene (riga)
mat_scaled <- t(scale(t(mat)))
top_genes <- apply(mat_scaled, 2, function(x) {
  names(sort(x, decreasing = TRUE))[1:10]
})
top_genes
#markers nemmeno per sbaglio cmq

#e invece per gene?
top_genes <- apply(assay(avg_expr), 1, function(x) {
  names(sort(x, decreasing = TRUE))[1]
})
marker_genes_list_singleR_annot <- readRDS("~/marker_genes_list_singleR_annot.rds")
markers <- unlist(marker_genes_list_singleR_annot)
top_genes[names(top_genes) %in% markers]
#fantastico proprio


#nuclei ------------------------------------------------------------------------------------------------------------
nuclei_list_ann <- readRDS("~/nuclei_list_ann.rds")

#i nuclei li devo ancora normalizzare
nuclei_list_ann <- lapply(nuclei_list_ann, function(spe){
  spe <- logNormCounts(spe)
})

#mettiamo assieme gli sham e partiamo da là
nuclei_list_ann <- lapply(nuclei_list_ann, function(spe) {
  keep <- !is.na(spe$cell_type) & colSums(counts(spe)) > 0
  spe[,keep]
})
for (nm in names(nuclei_list_ann)) {
  spe <- nuclei_list_ann[[nm]]
  colData(spe)$sample_id <- rep(nm, ncol(spe))
  nuclei_list_ann[[nm]] <- spe
}
nuclei_list_ann$blocco1_sham$in_treatment <- rep(NA,dim(nuclei_list_ann$blocco1_sham)[2])
nuclei_list_ann$blocco3_sham$in_treatment <- rep(NA,dim(nuclei_list_ann$blocco3_sham)[2])
nuclei_list_ann$blocco6_sham$in_treatment <- rep(NA,dim(nuclei_list_ann$blocco6_sham)[2])
nuclei_list_ann$blocco4_c26$in_treatment <- rep(FALSE,dim(nuclei_list_ann$blocco4_c26)[2])
nuclei_list_ann$blocco6_c26$in_treatment <- rep(FALSE,dim(nuclei_list_ann$blocco6_c26)[2])
nuclei_list_ann$blocco5_c26STAT3$in_treatment <- rep(TRUE,dim(nuclei_list_ann$blocco5_c26STAT3)[2])
combined_spe <- do.call(cbind, nuclei_list_ann)
sce_sham <- combined_spe[,combined_spe$sampletype == "sham"]

#average expression for each gene in each celltype
library(scuttle)
ids <- paste(
  sce_sham$sample_id,
  sce_sham$cell_type,
  sep = "_"
)
avg_expr <- aggregateAcrossCells(
  sce_sham,
  ids = ids,
  use.assay.type = "logcounts"
)

#number of detected genes
library(Matrix)
counts_mat <- assay(sce_sham, "counts") > 0
genes_per_cell <- colSums(counts_mat)
sce_sham$detected_genes <- genes_per_cell
tapply(sce_sham$detected_genes,
       sce_sham$cell_type,
       mean)
#for example
sce_sham[,sce_sham$cell_type == "MuSC"]$detected_genes

# most expressed genes for each cell type
# top_genes <- apply(assay(avg_expr), 2, function(x) {
#   names(sort(x, decreasing = TRUE))[1:10]
# })
# top_genes

#but maybe if I divided for libsize
mat <- assay(avg_expr)
# scaling per gene (riga)
mat_scaled <- t(scale(t(mat)))
top_genes <- apply(mat_scaled, 2, function(x) {
  names(sort(x, decreasing = TRUE))[1:10]
})
top_genes
#markers nemmeno per sbaglio cmq

#e invece per gene?
top_genes <- apply(assay(avg_expr), 1, function(x) {
  names(sort(x, decreasing = TRUE))[1]
})
marker_genes_list_singleR_annot <- readRDS("~/marker_genes_list_singleR_annot.rds")
markers <- unlist(marker_genes_list_singleR_annot)
top_genes[names(top_genes) %in% markers]
#fantastico proprio



















