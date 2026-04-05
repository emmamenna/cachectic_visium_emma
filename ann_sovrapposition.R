library(ggplot2)
library(ggnewscale)
library(SpatialExperiment)
library(scales)
library(viridis)
# 
# 
# #TRIM/MYH4 SOVRAPPOSITION WITH 16BINS AND NUCLEI
# nuclei_list <- readRDS("~/nuclei_list_ann.rds")
# c26_list <- nuclei_list[names(nuclei_list) %in% c("blocco4_c26","blocco6_c26")]
# 
# lapply(names(c26_list), function(name) {
#   
#   spe <- c26_list[[name]]
#   spe <- spe[,spe$to_discard == "FALSE"]
#   
#   trim_counts <- as.numeric(counts(spe["Trim63", ]))
#   ttn_counts <- as.numeric(counts(spe["Ttn", ]))
#   myh4_counts <- as.numeric(counts(spe["Myh4", ]))
#   myh1_counts <- as.numeric(counts(spe["Myh1", ]))
#   myh7_counts <- as.numeric(counts(spe["Myh7", ]))
#   
#   spe$trim <- (trim_counts > 0) & (ttn_counts > 0) & (myh4_counts == 0) &
#     (myh1_counts == 0) & (myh7_counts == 0)
#   
#   coords <- spatialCoords(spe)
#   
#   df <- data.frame(
#     x    = coords[, 2],
#     y    = coords[, 1],
#     trim = spe$trim
#   )
#   
#   df$trim_plot <- ifelse(df$trim, TRUE, NA)
# 
#   col_trim_true <- "orange"
#   
#   p <- ggplot(df, aes(x = x, y = y)) +
#     geom_point(
#       data = subset(df, trim_plot == TRUE),
#       aes(x = x, y = y),
#       shape = 16,
#       size  = 2,
#       color = col_trim_true
#     ) +
#     
#     coord_equal() +
#     scale_y_reverse() +
#     theme_void() +
#     theme(
#       legend.position = "right"
#     )
#   ggtitle(name)
#   p
# })
# 
# ##########################
# bin16_list <- readRDS("~/subset_list_updated.rds")
# c26_bin_list <- bin16_list[names(bin16_list) %in% c("c26_b4","c26_b6")]
# names(c26_bin_list) <- c("blocco4_c26","blocco6_c26")
# 
# lapply(names(c26_bin_list), function(name) {
#   
#   spe <- c26_bin_list[[name]]
#   myh4_counts <- as.numeric(counts(spe["Myh4", ]))
#   spe$myh4 <- (myh4_counts > 0)
#   
#   coords <- spatialCoords(spe)
#   
#   df <- data.frame(
#     x    = coords[, 1],
#     y    = coords[, 2],
#     trim = spe$myh4
#   )
#   
#   df$myh4 <- myh4_counts > 0
#   
#   p <- ggplot(df, aes(x = x, y = y)) +
#     geom_point(
#       data = df[df$myh4, ],
#       shape = 16,
#       size  = 2,
#       color = "lightblue"
#     ) +
#     coord_equal() +
#     scale_y_reverse() +
#     theme_void() +
#     ggtitle(name)
# 
# })
# 
# #UNICO GRAFICO -------------------------------------------------------------------
# library(ggplot2)
# library(patchwork)  # install.packages("patchwork") se non ce l'hai
# 
# lapply(names(c26_list), function(name) {
#   
#   ## --- TRIM plot ---
#   spe1 <- c26_list[[name]]
#   spe1 <- spe1[, spe1$to_discard == "FALSE"]
#   
#   trim_counts <- as.numeric(counts(spe1["Trim63", ]))
#   ttn_counts  <- as.numeric(counts(spe1["Ttn", ]))
#   myh4_counts1 <- as.numeric(counts(spe1["Myh4", ]))
#   myh1_counts <- as.numeric(counts(spe1["Myh1", ]))
#   myh7_counts <- as.numeric(counts(spe1["Myh7", ]))
#   
#   trim_flag <- (trim_counts > 0) & (ttn_counts > 0) &
#     (myh4_counts1 == 0) & (myh1_counts == 0) & (myh7_counts == 0)
#   
#   coords1 <- spatialCoords(spe1)
#   df1 <- data.frame(
#     x = coords1[, 2],
#     y = coords1[, 1],
#     trim = trim_flag
#   )
#   
#   p1 <- ggplot(df1, aes(x = x, y = y)) +
#     geom_point(data = df1[df1$trim, ], color = "orange", size = 2) +
#     coord_equal() +
#     scale_y_reverse() +
#     theme_void() +
#     ggtitle(paste0(name, " - TRIM"))
#   
#   ## --- MYH4 plot ---
#   spe2 <- c26_bin_list[[name]]
#   myh4_counts2 <- as.numeric(counts(spe2["Myh4", ]))
#   
#   coords2 <- spatialCoords(spe2)
#   df2 <- data.frame(
#     x = coords2[, 1],
#     y = coords2[, 2],
#     myh4 = myh4_counts2 > 0
#   )
#   
#   p2 <- ggplot(df2, aes(x = x, y = y)) +
#     geom_point(data = df2[df2$myh4, ], color = "lightblue", size = 2) +
#     coord_equal() +
#     scale_y_reverse() +
#     theme_void() +
#     ggtitle(paste0(name, " - MYH4"))
#   
#   ## --- Combina verticalmente ---
#   p_combined <- p1 / p2   # patchwork usa / per vertical, | per orizzontale
#   p_combined
# })


#Annotation based trim-2b sovrapposition ----------------------------------------------------------

nuclei_list <- readRDS("~/nuclei_list_ann.rds")
c26_list <- nuclei_list[names(nuclei_list) %in% c("blocco4_c26","blocco6_c26")]

default_color <- "white"
your_color_vector <- setNames(
  rep(default_color, 16),
  c("Endothelial", "FAPs", "Immune_Cells", "MuSC", "Myonuclei_IIx",
    "Myonuclei_IIx_IIa", "Myonuclei_IIx_IIb", "Myonuclei_MTJ",
    "Myonuclei_NMJ", "Nervous_System", "Pericyte",
    "Smooth_Muscular", "Tenocyte",
    "Myonuclei_IIb", "Myonuclei_Trim63")
)
your_color_vector["Myonuclei_Trim63"] <- "orange"

plot_list <- lapply(names(c26_list), function(nm) {
  spe <- c26_list[[nm]]
  spe <- spe[,spe$to_discard == "FALSE"]
  spe$in_tissue <- rep(TRUE,dim(spe)[2])
  plotCoords(spe, annotate = "cell_type", point_size = 0.9,
             x_coord = "y_coord", y_coord = "x_coord") +
    scale_color_manual(values = your_color_vector, na.value = default_color) +
    ggtitle(nm) + 
    theme(
      legend.key.width  = unit(0.5, "lines"),
      legend.key.height = unit(1, "lines"),
      plot.title = element_text(hjust = 0.5)
    )
})

bin16_list <- readRDS("~/subset_list_updated.rds")
c26_bin_list <- bin16_list[names(bin16_list) %in% c("c26_b4","c26_b6")]
names(c26_bin_list) <- c("blocco4_c26","blocco6_c26")

default_color <- "white"
your_color_vector <- setNames(
  rep(default_color, 16),
  c("Endothelial", "FAPs", "Immune_Cells", "MuSC", "Myonuclei_IIx",
    "Myonuclei_IIx_IIa", "Myonuclei_IIx_IIb", "Myonuclei_MTJ",
    "Myonuclei_NMJ", "Nervous_System", "Pericyte",
    "Smooth_Muscular", "Tenocyte",
    "Myonuclei_IIb", "Myonuclei_Trim63")
)
your_color_vector["Myonuclei_IIb"] <- "lightblue"
plot_list <- lapply(names(c26_bin_list), function(nm) {
  spe <- c26_bin_list[[nm]]
  spe$in_tissue <- rep(TRUE,dim(spe)[2])
  plotCoords(spe, annotate = "new_cell_type", point_size = 0.9) +
    scale_color_manual(values = your_color_vector, na.value = default_color) +
    ggtitle(nm) + 
    theme(
      legend.key.width  = unit(0.5, "lines"),
      legend.key.height = unit(1, "lines"),
      plot.title = element_text(hjust = 0.5)
    )
})

#GRAFICO UNICO ------------------------------------------------------------------------------

range(spatialCoords(c26_list$blocco4_c26)[,1])
range(spatialCoords(c26_list$blocco4_c26)[,2])

range(spatialCoords(c26_bin_list$blocco4_c26)[,2])
range(spatialCoords(c26_bin_list$blocco4_c26)[,1])

library(ggplot2)
library(SingleCellExperiment)

plot_list <- lapply(names(c26_list), function(nm) {
  
  spe_nuclei <- c26_list[[nm]][, c26_list[[nm]]$to_discard == "FALSE"]
  spe_bin    <- c26_bin_list[[nm]]
  
  # Estrai coordinate e annotazioni
  df_bin <- data.frame(
    x = scale(spatialCoords(spe_bin)[,1]),
    y = scale(spatialCoords(spe_bin)[,2]),
    new_cell_type = colData(spe_bin)$new_cell_type
  )
  print(range(df_bin$y))
  df_bin$color <- ifelse(df_bin$new_cell_type == "Myonuclei_IIb", "lightblue", "white")
  
  df_nuclei <- data.frame(
    x = scale(spatialCoords(spe_nuclei)[,2]),
    y = scale(spatialCoords(spe_nuclei)[,1]),
    cell_type = colData(spe_nuclei)$cell_type
  )
  print(range(df_nuclei$y))
  df_nuclei$color <- ifelse(df_nuclei$cell_type == "Myonuclei_Trim63", "orange", "white")
  
  # Combina i due data.frame in uno solo
  df_plot <- rbind(
    data.frame(x = df_bin$x, y = df_bin$y, color = df_bin$color),
    data.frame(x = df_nuclei$x, y = df_nuclei$y, color = df_nuclei$color)
  )
  
  # Plot unico con un solo mapping colore
  ggplot(df_plot, aes(x = x, y = y, color = color)) +
    geom_point(size = 0.9, alpha = 0.8) +
    scale_color_identity() +  # usa direttamente i colori già definiti
    ggtitle(nm) +
    theme(
      legend.position = "none",
      plot.title = element_text(hjust = 0.5)
    )
})
plot_list

##########################################################
plot_list <- lapply(names(c26_list), function(nm) {
  
  spe_nuclei <- c26_list[[nm]][, c26_list[[nm]]$to_discard == "FALSE"]
  spe_bin    <- c26_bin_list[[nm]]
  # Define cell types to keep
  bin_types <- c("Endothelial", "FAPs", "Immune_Cells", "MuSC", "Myonuclei_IIx",
                 "Myonuclei_IIx_IIa", "Myonuclei_IIx_IIb", "Myonuclei_MTJ",
                 "Myonuclei_NMJ", "Nervous_System", "Pericyte",
                 "Smooth_Muscular", "Tenocyte",
                 "Myonuclei_IIb", "Myonuclei_Trim63")
  
  df_bin <- data.frame(
    x = scale(spatialCoords(spe_bin)[,1]),
    y = scale(spatialCoords(spe_bin)[,2]),
    new_cell_type = colData(spe_bin)$new_cell_type
  )
  
  # Keep only selected types
  df_bin <- df_bin[df_bin$new_cell_type %in% bin_types, ]
  
  # Assign colors
  df_bin$color <- ifelse(
    df_bin$new_cell_type == "Myonuclei_IIb",
    "lightblue",
    "lightgrey"
  )
  
  df_nuclei <- data.frame(
    x = scale(spatialCoords(spe_nuclei)[,2]),
    y = scale(spatialCoords(spe_nuclei)[,1]),
    cell_type = colData(spe_nuclei)$cell_type
  )
  df_nuclei <- df_nuclei[df_nuclei$cell_type == "Myonuclei_Trim63", ]
  df_nuclei$color <- "orange"
  
  df_plot <- rbind(
    data.frame(x = df_bin$x, y = df_bin$y, color = df_bin$color),
    data.frame(x = df_nuclei$x, y = df_nuclei$y, color = df_nuclei$color)
  )
  
  ggplot(df_plot, aes(x = x, y = y, color = color)) +
    geom_point(size = 0.9, alpha = 0.8) +
    scale_color_identity() +
    ggtitle(nm) +
    theme_void() +  # rimuove assi, griglie e sfondo
    theme(
      legend.position = "none",
      plot.title = element_text(hjust = 0.5),
      panel.background = element_blank(),
      panel.grid = element_blank()
    )
})

plot_list

library(purrr)
nuclei_list_ann <- readRDS("~/nuclei_list_ann.rds")
bin16_list_ann <- readRDS("~/bin16_list_ann.rds")

df_annot_nuclei <- map(nuclei_list_ann, function(spe_n){
  data.frame(
    nuclei_names = colnames(spe_n),
    nuclei_types = colData(spe_n)$cell_type
  )
})

df_annot_bin <- map(bin16_list_ann, function(spe_b){
  data.frame(
    bin_names = colnames(spe_b),
    bin_types = colData(spe_b)$new_cell_type
  )
})

library(dplyr)
library(purrr)
library(arrow)

df_combined <- imap_dfr(df_annot_bin, ~ .x %>%
                          mutate(sample = .y))
write_parquet(df_combined, "bin_plot_info.parquet")
df_combined_nuclei <- imap_dfr(df_annot_nuclei, ~ .x %>%
                          mutate(sample = .y))
write_parquet(df_combined_nuclei, "nuclei_plot_info.parquet")




