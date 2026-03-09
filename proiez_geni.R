load("~/subset_list.RData")

library(ggplot2)
library(SpatialExperiment)
library(patchwork)

cac_list <- subset_list[names(subset_list) %in% c("c26_b2","c26_b4","c26_b6")]
geni <- c("Bax","Casp9","Casp8","Casp3","Apaf1","Trp53","Pax7","Myod1","Myf5","Myh3")
#"Bel2"??

library(ggplot2)
library(SpatialExperiment) # assuming your objects are from this package

# cac_list: list of SpatialExperiment objects
# geni: vector of gene names to plot

for (spe_name in names(cac_list)) {
  spe <- cac_list[[spe_name]] # Use named access
  
  # Get main counts matrix
  expr_matrix <- counts(spe)[geni, , drop = FALSE]
  
  coords <- as.data.frame(spatialCoords(spe))
  colnames(coords) <- c("x_coord", "y_coord")
  
  for (i in seq_along(geni)) {
    marker <- geni[i]
    expr <- as.numeric(expr_matrix[marker, ])
    
    coords$expr <- expr
    
    p <- ggplot(coords, aes(x = x_coord, y = y_coord, color = expr)) +
      geom_point(size = 0.5) +
      coord_fixed() +
      scale_y_reverse() +
      ggtitle(paste0(spe_name, ": ", marker)) +
      theme_void() +
      theme(plot.title = element_text(size = 10)) +
      scale_color_gradientn(colors = rev(hcl.colors(9, "Rocket")))
    
    # Save the plot
    ggsave(filename = paste0(spe_name, "_", marker, ".png"), plot = p, width = 4, height = 4, dpi = 300)
  }
}


