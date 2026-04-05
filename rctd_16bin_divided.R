library(SingleCellExperiment)
library(SpatialExperiment)
library(SummarizedExperiment)
library(spacexr)
library(dplyr)
library(ggplot2)
library(Matrix)
library(S4Vectors)

sce_ref <- readRDS("~/multiome_sce.rds") #reference dataset from multiome
subset_list <- readRDS("~/subset_list.rds")

color_vector <- c(
  "Endothelial" = "#0072B2",       # blu vivo
  "FAPs" = "#009E73",              # verde vivo
  "Immune_Cells" = "#D55E00",      # arancione-rosso
  "MuSC" = "#E69F00",              # arancione
  "Myonuclei_IIb" = "#CC79A7",     # magenta
  "Myonuclei_IIx" = "#56B4E9",     # azzurro chiaro
  "Myonuclei_IIx_IIa" = "#F0E442", # giallo
  "Myonuclei_IIx_IIb" = "#9B59B6", # viola chiaro
  "Myonuclei_MTJ" = "#0099B4",     # turchese
  "Myonuclei_NMJ" = "#DDAA33",     # senape
  "Myonuclei_Trim63" = "#CC3311",  # rosso scuro
  "Nervous_System" = "#44AA99",    # verde acqua
  "Pericyte" = "#AA4499",          # viola/rosa
  "Smooth_Muscular" = "#332288",   # blu-viola scuro
  "Tenocyte" = "#A3E635"           # verde lime
)

score_cols_sham <- c("Endothelial", "FAPs", "Immune_Cells", "MuSC", 
                "Myonuclei_IIb", "Myonuclei_IIx", "Myonuclei_IIx_IIa", 
                "Myonuclei_IIx_IIb", "Myonuclei_MTJ", "Myonuclei_NMJ",
                "Nervous_System", "Pericyte", 
                "Smooth_Muscular", "Tenocyte")

score_cols_cac <- c("Endothelial", "FAPs", "Immune_Cells", "MuSC", 
                "Myonuclei_IIb", "Myonuclei_IIx", "Myonuclei_IIx_IIa", 
                "Myonuclei_IIx_IIb", "Myonuclei_MTJ", "Myonuclei_NMJ", "Myonuclei_Trim63",
                "Nervous_System", "Pericyte", 
                "Smooth_Muscular", "Tenocyte")

#FUNCTIONS ------------------------------------------------------------------------------

run_my_rctd <- function(spe, sce_ref, out_prefix) {
  set.seed(9237)
  rownames(spe) <- rowData(spe)$Symbol
  rownames(spe) <- make.unique(rownames(spe))
  
  #data preparation
  rctd_data <- createRctd(
    spatial_experiment = spe,
    reference_experiment = sce_ref,
    cell_type_col = "subclusters",
    ref_n_cells_min = 15
  )
  
  # For multi mode:
  # cell_type_list: List of cell types per pixel
  # conf_list: List of whether cell type predictions are confident
  # Additional metrics like min_score (likelihood)
  
  #deconvolution
  results <- runRctd(rctd_data, rctd_mode = "multi")
  
  #weights extraction
  ws <- assay(results, "weights")
  ws <- data.frame(t(as.matrix(ws)))
  colData(spe)[names(ws)] <- ws[colnames(spe), ] 
  
  #Cluster labelling
  table(rowSums(ws == 1))
  ids <- names(ws)[apply(ws, 1, which.max)]
  ids <- gsub("\\.([A-z])", " \\1", ids)
  idx <- match(colnames(spe), rownames(ws))
  spe$cell_type2 <- factor(ids[idx]) 
  
  spe
}

fix_spe_celltypes <- function(spe, score_columns, cell_type_col = "cell_type") {
  cd <- as.data.frame(colData(spe))
  score_mat <- as.matrix(cd[, score_columns])
  max_score <- apply(score_mat, 1, max)
  best_type <- score_columns[max.col(score_mat, ties.method = "first")]
  orig_col_idx <- match(as.character(cd[[cell_type_col]]), score_columns)
  orig_score <- score_mat[cbind(1:nrow(score_mat), orig_col_idx)]
  new_cell_type <- as.character(cd[[cell_type_col]])
  reassign_mask <- which(orig_score == 0 & max_score > 0)
  new_cell_type[reassign_mask] <- best_type[reassign_mask]
  na_mask <- which(max_score == 0)
  new_cell_type[na_mask] <- NA
  spe$new_cell_type <- new_cell_type
  return(spe)
}

#ONLY SHAM ON SHAM ---------------------------------------------------------------------------------
sham_ref <- sce_ref[,sce_ref$condition=="SHAM"]
sham_ref <- sham_ref[,sham_ref$subclusters != "Myonuclei_Trim63"]

#block1
shamb1 <- subset_list$sham_b1
shamb1 <- run_my_rctd(shamb1, sham_ref, "sham_b1")
shamb1 <- fix_spe_celltypes(shamb1,score_cols_sham,cell_type_col = "cell_type2")
shamb1 <- shamb1[,!is.na(shamb1$new_cell_type)]

# plotCoords(shamb1,annotate="new_cell_type",point_size = 0.5) +
#   scale_color_manual(values = color_vector) +
#   theme(
#     legend.key.width = unit(0.5, "lines"),
#     legend.key.height = unit(1, "lines")
#   )
# round(prop.table(table(shamb1$new_cell_type)),3)

#block3
shamb3 <- subset_list$sham_b3
shamb3 <- run_my_rctd(shamb3, sham_ref, "sham_b3")
shamb3 <- fix_spe_celltypes(shamb3,score_cols_sham, cell_type_col = "cell_type2")
shamb3 <- shamb3[,!is.na(shamb3$new_cell_type)]

# plotCoords(shamb3,annotate="new_cell_type",point_size = 0.9) +
#   scale_color_manual(values = color_vector) +
#   theme(
#     legend.key.width = unit(0.5, "lines"),
#     legend.key.height = unit(1, "lines")
#   )
# round(prop.table(table(shamb3$new_cell_type)),3)

#block6
shamb6 <- subset_list$sham_b6
shamb6 <- run_my_rctd(shamb6, sham_ref, "sham_b6")
shamb6 <- fix_spe_celltypes(shamb6,score_cols_sham,cell_type_col = "cell_type2")
shamb6 <- shamb6[,!is.na(shamb6$new_cell_type)]

# plotCoords(shamb6,annotate="new_cell_type",point_size = 0.9) +
#   scale_color_manual(values = color_vector) +
#   theme(
#     legend.key.width = unit(0.5, "lines"),
#     legend.key.height = unit(1, "lines")
#   )
# round(prop.table(table(shamb6$new_cell_type)),3)

subset_list_updated <- readRDS("~/subset_list_updated.rds")
dim(shamb1)
dim(subset_list_updated$sham_b1)
dim(shamb3)
dim(subset_list_updated$sham_b3)
dim(shamb6)
dim(subset_list_updated$sham_b6)

#SHAM + trim on C26  ---------------------------------------------------------------------------------
cac_ref_sham <- sce_ref[,sce_ref$condition=="SHAM" |
                     (sce_ref$condition=="CAC" & sce_ref$subclusters == "Myonuclei_Trim63")]

#block4
c26_b4_sham <- subset_list$c26_b4
c26_b4_sham <- run_my_rctd(c26_b4_sham, cac_ref_sham, "c26_b4")
c26_b4_sham <- fix_spe_celltypes(c26_b4_sham,score_cols_cac,cell_type_col = "cell_type2")
c26_b4_sham <- c26_b4_sham[,!is.na(c26_b4_sham$new_cell_type)]

# plotCoords(c26_b4,annotate="new_cell_type",point_size = 0.9) +
#   scale_color_manual(values = color_vector) +
#   theme(
#     legend.key.width = unit(0.5, "lines"),
#     legend.key.height = unit(1, "lines")
#   )
# round(prop.table(table(c26_b4$new_cell_type)),3)

#block2
c26_b2_sham <- subset_list$c26_b2
c26_b2_sham <- run_my_rctd(c26_b2_sham, cac_ref_sham, "c26_b2")
c26_b2_sham <- fix_spe_celltypes(c26_b2_sham,score_cols_cac,cell_type_col = "cell_type2")
c26_b2_sham <- c26_b2_sham[,!is.na(c26_b2_sham$new_cell_type)]

# plotCoords(c26_b2_sham,annotate="new_cell_type",point_size = 0.9) +
#   scale_color_manual(values = color_vector) +
#   theme(
#     legend.key.width = unit(0.5, "lines"),
#     legend.key.height = unit(1, "lines")
#   )
# round(prop.table(table(c26_b2_sham$new_cell_type)),3)

#block6
c26_b6_sham <- subset_list$c26_b6
c26_b6_sham <- run_my_rctd(c26_b6_sham, cac_ref_sham, "c26_b6")
c26_b6_sham <- fix_spe_celltypes(c26_b6_sham,score_cols_cac,cell_type_col = "cell_type2")
c26_b6_sham <- c26_b6_sham[,!is.na(c26_b6_sham$new_cell_type)]

# plotCoords(c26_b6_sham, annotate="new_cell_type",point_size = 0.9) +
#   scale_color_manual(values = color_vector) +
#   theme(
#     legend.key.width = unit(0.5, "lines"),
#     legend.key.height = unit(1, "lines")
#   )
# round(prop.table(table(c26_b6_sham$new_cell_type)),3)


#ONLY CAC on C26  ---------------------------------------------------------------------------------
cac_ref <- sce_ref[,sce_ref$condition=="CAC"]

#block4
c26_b4 <- subset_list$c26_b4
#table(cac_ref$subclusters)
#diversi cell types hanno numerosità troppo bassa così
#ho cambiato il parametro della funzione
c26_b4 <- run_my_rctd(c26_b4, cac_ref, "c26_b4")
c26_b4 <- fix_spe_celltypes(c26_b4,score_cols_cac,cell_type_col = "cell_type2")
c26_b4 <- c26_b4[,!is.na(c26_b4$new_cell_type)]

# plotCoords(c26_b4,annotate="new_cell_type",point_size = 0.9) +
#   scale_color_manual(values = color_vector) +
#   theme(
#     legend.key.width = unit(0.5, "lines"),
#     legend.key.height = unit(1, "lines")
#   )
# round(prop.table(table(c26_b4$new_cell_type)),3)

#block2
c26_b2 <- subset_list$c26_b2
c26_b2 <- run_my_rctd(c26_b2, cac_ref, "c26_b2")
table(c26_b2$cell_type2)
c26_b2 <- fix_spe_celltypes(c26_b2,score_cols_cac,cell_type_col = "cell_type2")
keep <- !is.na(c26_b2$new_cell_type)
c26_b2 <- c26_b2[,keep]

# plotCoords(c26_b2,annotate="new_cell_type",point_size = 0.9) +
#   scale_color_manual(values = color_vector) +
#   theme(
#     legend.key.width = unit(0.5, "lines"),
#     legend.key.height = unit(1, "lines")
#   )
# round(prop.table(table(c26_b2$new_cell_type)),3)

#block 6
c26_b6 <- subset_list$c26_b6
c26_b6 <- run_my_rctd(c26_b6, cac_ref, "c26_b6")
c26_b6 <- fix_spe_celltypes(c26_b6,score_cols,cell_type_col = "cell_type2")
keep <- !is.na(c26_b6$new_cell_type)
c26_b6 <- c26_b6[,keep]

# plotCoords(c26_b6,annotate="new_cell_type",point_size = 0.9) +
#   scale_color_manual(values = color_vector) +
#   theme(
#     legend.key.width = unit(0.5, "lines"),
#     legend.key.height = unit(1, "lines")
#   )
# round(prop.table(table(c26_b6$new_cell_type)),3)

#OVERALL PROPORTIONS ------------------------------------------------------------------------------
sham_list <- list(shamb1,shamb3,shamb6)

prop_df_condition <- lapply(sham_list, function(spe) {
  data.frame(
    cell_type = spe$new_cell_type
  )
}) |> bind_rows() |> filter(!is.na(cell_type))

prop_df_condition <- prop_df_condition %>%
  dplyr::count(cell_type) %>%
  dplyr::mutate(proportion = n / sum(n))

ggplot(prop_df_condition, aes(x = "", y = proportion, fill = cell_type)) +
  geom_bar(stat = "identity", width = 0.5) +
  scale_fill_manual(values = color_vector) +
  theme_minimal() +
  labs(
    title = "Proporzioni dei tipi cellulari",
    x = NULL,
    y = "Proporzione",
    fill = "Tipo Cellulare"
  )

#CACHECTIC
c26_cac_list <- list(c26_b2_sham,c26_b4_sham,c26_b6_sham)

prop_df_condition <- lapply(c26_cac_list, function(spe) {
  data.frame(
    cell_type = spe$new_cell_type
  )
}) |> bind_rows() |> filter(!is.na(cell_type))

prop_df_condition <- prop_df_condition %>%
  dplyr::count(cell_type) %>%
  dplyr::mutate(proportion = n / sum(n))

ggplot(prop_df_condition, aes(x = "", y = proportion, fill = cell_type)) +
  geom_bar(stat = "identity", width = 0.5) +
  scale_fill_manual(values = color_vector) +
  theme_minimal() +
  labs(
    title = "Proporzioni dei tipi cellulari",
    x = NULL,
    y = "Proporzione",
    fill = "Tipo Cellulare"
  )

