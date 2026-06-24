subset_list <- readRDS("bin16_list_ann.rds")

library(SpatialExperiment)
library(ggplot2)
library(ggspavis)

# for (spe_name in names(subset_list)) {
#   spe <- subset_list[[spe_name]]
#   dummies <- model.matrix(~ spe$new_cell_type - 1)
#   
#   stat1 <- plotCoords(spe, annotate = "", point_size = 0.5) +
#     scale_color_manual(values = color_vector) +
#     theme(
#       legend.key.width = unit(0.5, "lines"),
#       legend.key.height = unit(1, "lines")
#     )
# }

subset_list <- lapply(subset_list, function(spe){
  spe <- spe[,!is.na(spe$new_cell_type)]
})

for (spe_name in names(subset_list)) {
  spe <- subset_list[[spe_name]]
  levels_list <- levels(as.factor(spe$new_cell_type))
  
  for (lev in levels_list) {
    spe$annotation_label <- ifelse(spe$new_cell_type == lev, "yes", "no")
    spe$annotation_label <- as.factor(spe$annotation_label)

    plot <- plotCoords(
      spe,
      annotate = "annotation_label", 
      point_size = 0.5
    ) +
      scale_color_manual(values = c("no" = "grey80", "yes" = "violet")) +
      theme(
        legend.key.width = unit(0.5, "lines"),
        legend.key.height = unit(1, "lines")
      ) +
      ggtitle(paste(lev,"cells in",spe_name))
    
    ggsave(paste0("fig_", spe_name, "_", lev, ".png"), plot)
  }
}

#per i nuclei
nuclei_list_ann <- readRDS("~/nuclei_list_ann.rds")
spe <- nuclei_list_ann$blocco2_c26
spe$trim <- ifelse(spe$cell_type == "Myonuclei_Trim63", "yes", "no")
spe$trim <- as.factor(spe$trim)
spe$in_tissue <- rep(TRUE,dim(spe)[2])
plot <- plotCoords(spe, annotate = "trim",  point_size = 0.5,
                   x_coord = "y_coord", y_coord = "x_coord") +
  scale_color_manual(values = c("no" = "grey80", "yes" = "red")) +
  theme(
    legend.key.width = unit(0.5, "lines"),
    legend.key.height = unit(1, "lines")
  ) 
ggsave("nuclei_trim.png", plot)



