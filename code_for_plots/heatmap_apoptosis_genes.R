library(SingleCellExperiment)
library(SummarizedExperiment)
library(scuttle)
library(ggplot2)
library(dplyr)
library(ComplexHeatmap)

subset_list <- readRDS("~/bin16_list_ann.rds")
#sample identifiers
for (nm in names(subset_list)) {
  spe <- subset_list[[nm]]
  colData(spe)$sample_id <- rep(nm, ncol(spe))
  subset_list[[nm]] <- spe
}
subset_list$c26SMAD23_b2$to_discard <- NULL
subset_list$c26murf1_b2$to_discard <- NULL
combined_spe <- do.call(cbind, subset_list)
table(combined_spe$tissue_type)

#1. trim c26 vs myo overall sham ---------------------------------------------------------------
genes_apop_over <- c("Ctsl", "Tuba1c", "Gadd45a","Bcl2l1", "Ern1", "Ctsd", "Nol3", "Daxx", "Rela", "Akt1",
                      "Pik3r1", "Ctss", "Ctsb", "Parp4", "Gadd45b", "Fas", "Bak1", "Map3k14", "Ctsz", 
                      "Ripk1", "Tnfrsf1a", "Pdpk1", "Eif2s1", "Atf4", "Diablo", "Mcl1")
spe_sub <- combined_spe[, (combined_spe$tissue_type == "c26" & 
                           !is.na(combined_spe$new_cell_type) &
                           combined_spe$new_cell_type == "Myonuclei_Trim63") |
                          (combined_spe$tissue_type == "sham" & 
                             !is.na(combined_spe$new_cell_type) &
                             combined_spe$new_cell_type %in% c(
                               "Myonuclei_IIb", "Myonuclei_IIx", "Myonuclei_IIx_IIa",
                               "Myonuclei_IIx_IIb", "Myonuclei_MTJ",
                               "Myonuclei_NMJ", "Myonuclei_Trim63"
                             ))
]
mat <- logcounts(spe_sub)[genes_apop_over, ]
group <- ifelse(
  spe_sub$tissue_type == "c26" & 
    spe_sub$new_cell_type == "Myonuclei_Trim63",
  "c26_Trim63",
  "sham_myo_overall"
)
mat_avg <- sapply(unique(group), function(g) {
  rowMeans(mat[, group == g, drop = FALSE])
})

ht <- Heatmap(
  mat_avg,
  name = "expr",
  
  cluster_rows = TRUE,
  cluster_columns = FALSE,
  
  column_title = "trim_c26 vs myo_sham",
  
  row_names_gp = grid::gpar(fontsize = 8)
)
pdf("heatmap_overall.pdf",
    width = 6,
    height = 10,
    useDingbats = FALSE)
draw(ht)
dev.off()

#2. trim c26 vs myo 2b sham ---------------------------------------------------------------
genes_apop_2b <- c( "Ctsl", "Tuba1c", "Ctss", "Bcl2l1", "Ern1", "Gadd45a", "Akt1", "Pik3r1",
                    "Rela", "Fas", "Parp4", "Ctsd", "Nol3", "Daxx", "Map3k14", "Tnfrsf1a", "Ripk1",
                    "Ctsb", "Ctsf", "Ctsz", "Mcl1", "Ctsh", "Sptan1", "Diablo", "Gadd45b", "Eif2s1",
                    "Ikbkg", "Fadd")
  
spe_sub <- combined_spe[, (combined_spe$tissue_type == "c26" & 
                             !is.na(combined_spe$new_cell_type) &
                             combined_spe$new_cell_type == "Myonuclei_Trim63") |
                          (combined_spe$tissue_type == "sham" & 
                             !is.na(combined_spe$new_cell_type) &
                             combined_spe$new_cell_type == "Myonuclei_IIb")
]
mat <- logcounts(spe_sub)[genes_apop_2b, ]
group <- ifelse(
  spe_sub$tissue_type == "c26" & 
    spe_sub$new_cell_type == "Myonuclei_Trim63",
  "c26_Trim63",
  "sham_myo_2b"
)
mat_avg <- sapply(unique(group), function(g) {
  rowMeans(mat[, group == g, drop = FALSE])
})

ht <- Heatmap(
  mat_avg,
  name = "expr",
  
  cluster_rows = TRUE,
  cluster_columns = FALSE,
  
  column_title = "trim_c26 vs 2b_sham",
  
  row_names_gp = grid::gpar(fontsize = 8)
)
pdf("heatmap_2b.pdf",
    width = 6,
    height = 10,
    useDingbats = FALSE)
draw(ht)
dev.off()

