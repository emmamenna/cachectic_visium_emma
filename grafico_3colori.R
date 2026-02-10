# load("~/subset_list.RData")
# c26_list <- subset_list[names(subset_list) %in% c("c26_b2","c26_b4","c26_b6")]

#con gli 8bin
library(ggplot2)
library(ggspavis)
library(SpatialExperiment)
load("~/spe_list_8bin.RData")
names(spe_list)
c26sham_list <- spe_list[names(spe_list) %in% c("blocco1_sham","blocco2_c26","blocco3_sham",
                                                "blocco4_c26","blocco6_c26","blocco6_sham")]

#FAPS ---------------------------------------------------------------------------------
out_dir <- "faps_plots"
dir.create(out_dir, showWarnings = FALSE)

default_color <- "lightgrey" #+chiaro
your_color_vector <- setNames(
  rep(default_color, 2),
  c("TRUE","FALSE")
)
your_color_vector["TRUE"]  <- "blue"
genes_faps <- c("Pdgfra", "Ly6a", "Dcn")

lapply(names(c26sham_list), function(name) {
  
  spe <- c26sham_list[[name]]
  spe$faps <- colSums(counts(spe[genes_faps, ])) > 0
  
  p <- plotCoords(spe, annotate = "faps", point_size = 0.7) +
    scale_color_manual(values = your_color_vector, na.value = default_color) +
    theme(
      legend.key.width  = unit(0.5, "lines"),
      legend.key.height = unit(1, "lines"),
      plot.title = element_text(hjust = 0.5)
    ) +
    ggtitle(name)
  
  ggsave(
    filename = file.path(out_dir, paste0(name, ".png")),
    plot = p,
    width = 6,
    height = 5,
    dpi = 300
  )
})

#TRIM --------------------------------------------------------------------------------
out_dir <- "plots_trim"
dir.create(out_dir, showWarnings = FALSE)
lapply(names(c26sham_list), function(name) {
  
  spe <- c26sham_list[[name]]
  trim_counts <- counts(spe)["Trim63", ]
  spe$trim <- trim_counts
  p <- plotCoords(spe, annotate = "Trim63", point_size = 0.7) + ggtitle(name)
  
  ggsave(
    filename = file.path(out_dir, paste0(name, ".png")),
    plot = p,
    width = 6,
    height = 5,
    dpi = 300
  )
})

#MYH4 --------------------------------------------------------------------------------
out_dir <- "plots_myh4"
dir.create(out_dir, showWarnings = FALSE)

lapply(names(c26sham_list), function(name) {
  
  spe <- c26sham_list[[name]]
  myh4_counts <- counts(spe)["Myh4", ]
  spe$myh4 <- myh4_counts
  p <- plotCoords(spe, annotate = "myh4", point_size = 0.7) + ggtitle(name)
  
  ggsave(
    filename = file.path(out_dir, paste0(name, ".png")),
    plot = p,
    width = 6,
    height = 5,
    dpi = 300
  )
})

#PROVO A SOVRAPPORLI -------------------------------------------------------------------
library(ggplot2)
library(ggnewscale)
library(SpatialExperiment)
library(scales)
library(viridis)

out_dir <- "plots_combined"
dir.create(out_dir, showWarnings = FALSE)

genes_faps <- c("Pdgfra", "Ly6a", "Dcn")

lapply(names(c26sham_list), function(name) {
  
  spe <- c26sham_list[[name]]
  
  ## =========================
  ## 1. Annotazioni
  ## =========================
  spe$faps <- colSums(counts(spe[genes_faps, ])) > 0
  spe$trim <- counts(spe)["Trim63", ]
  spe$myh4 <- counts(spe)["Myh4", ]
  
  spe$trim[spe$trim == 0] <- NA
  spe$myh4[spe$myh4 == 0] <- NA
  
  ## =========================
  ## 2. Coordinate spaziali
  ## =========================
  coords <- spatialCoords(spe)
  
  df <- data.frame(
    x    = coords[, 1],
    y    = coords[, 2],
    faps = spe$faps,
    trim = spe$trim,
    myh4 = spe$myh4
  )
  
  ## =========================
  ## 3. Palette colori
  ## =========================
  col_faps_true  <- "orange"        
  col_faps_false <- "white"
  
  ## =========================
  ## 4. Plot
  ## =========================
  p <- ggplot(df, aes(x = x, y = y)) +
    
    ## ---- FAPs (sfondo)
    geom_point(
      aes(color = faps),
      size = 0.1
    ) +
    scale_color_manual(
      values = c("TRUE" = col_faps_true,
                 "FALSE" = col_faps_false),
      name = "FAPs",
      na.value = col_faps_false
    ) +
    
    ggnewscale::new_scale_color() +
    
    ## ---- Trim63 (verde-giallo)
    geom_point(
      aes(color = trim),
      size = 0.1
    ) +
    scale_color_gradientn(
      colors = viridis(100, option = "C"),
      name = "Trim63"
    ) +
    
    ggnewscale::new_scale_color() +
    
    ## ---- Myh4 (viola-arancio)
    geom_point(
      aes(color = myh4),
      size = 0.1
    ) +
    scale_color_gradientn(
      colors = viridis(100, option = "D"),
      name = "Myh4"
    ) +
    
    coord_equal() +
    scale_y_reverse() +
    theme_void() +
    ggtitle(name)
  
  ## =========================
  ## 5. Salvataggio
  ## =========================
  ggsave(
    filename = file.path(out_dir, paste0(name, "_combined.png")),
    plot = p,
    width = 6,
    height = 5,
    dpi = 300
  )
})


#SOLO FAPS E TRIM PER I C26 --------------------------------------------------------------
library(ggplot2)
library(ggnewscale)
library(SpatialExperiment)
library(scales)
library(viridis)

#c26_list <- spe_list[names(spe_list) %in% c("blocco2_c26","blocco4_c26","blocco6_c26")]

out_dir <- "plots_doppi"
dir.create(out_dir, showWarnings = FALSE)

genes_faps <- c("Pdgfra", "Ly6a", "Dcn")

lapply(names(c26_list), function(name) {
  
  spe <- c26_list[[name]]
  
  spe$faps <- colSums(counts(spe[genes_faps, ])) > 0
  spe$trim <- counts(spe)["Trim63", ]
  spe$myh4 <- counts(spe)["Myh4", ]
  
  spe$trim[spe$trim == 0] <- NA
  spe$myh4[spe$myh4 == 0] <- NA
  
  coords <- spatialCoords(spe)
  
  df <- data.frame(
    x    = coords[, 1],
    y    = coords[, 2],
    faps = spe$faps,
    trim = spe$trim,
    myh4 = spe$myh4
  )
  
  df$faps_plot <- ifelse(df$faps, TRUE, NA)
  df$trim_plot <- log1p(df$trim)
  df$myh4_plot <- log1p(df$myh4)
  
  col_faps_true <- "green"  # arancione caldo
  
  p <- ggplot(df, aes(x = x, y = y)) +
  geom_point(
    data = subset(df, !is.na(trim_plot)),
    aes(color = trim_plot),
    size = 0.2
  ) +
    scale_color_gradientn(
      colors = viridis(100, option = "C"),
      name = "Trim63 (log1p)"
    ) +
    
    ggnewscale::new_scale_color() +
  
    geom_point(
      data = subset(df, faps_plot == TRUE),
      aes(x = x, y = y),
      shape = 16,         # cerchio pieno
      size  = 0.4,
      color = "green"  # arancione pieno
    )+
    
    coord_equal() +
    scale_y_reverse() +
    theme_void() +
    ggtitle(name)
  
  ggsave(
    filename = file.path(out_dir, paste0(name, "_combined.png")),
    plot = p,
    width = 6,
    height = 5,
    dpi = 300
  )
})


#TRIM MENO MYH4 E FAPS ---------------------------------------------------------------
library(ggplot2)
library(ggnewscale)
library(SpatialExperiment)
library(scales)
library(viridis)

#c26_list <- spe_list[names(spe_list) %in% c("blocco2_c26","blocco4_c26","blocco6_c26")]

out_dir <- "plots_solotrim_faps"
dir.create(out_dir, showWarnings = FALSE)
  
  genes_faps <- c("Pdgfra", "Ly6a", "Dcn")
  
  lapply(names(c26_list), function(name) {
    
    spe <- c26_list[[name]]
    
    # =====================
    # 1. Annotazioni boolean
    # =====================
    spe$faps <- colSums(counts(spe[genes_faps, ])) > 0
    
    # Trim63 presente ma solo se Myh4 assente
    trim_counts <- as.numeric(counts(spe["Trim63", ]))
    myh4_counts <- as.numeric(counts(spe["Myh4", ]))
    
    spe$trim <- (trim_counts > 0) & (myh4_counts == 0)
    
    
    # =====================
    # 2. Coordinate
    # =====================
    coords <- spatialCoords(spe)
    
    df <- data.frame(
      x    = coords[, 1],
      y    = coords[, 2],
      faps = spe$faps,
      trim = spe$trim
    )
    
    # solo TRUE → NA per gli altri
    df$faps_plot <- ifelse(df$faps, TRUE, NA)
    df$trim_plot <- ifelse(df$trim, TRUE, NA)
    
    # =====================
    # 3. Palette
    # =====================
    col_faps_true <- "green"
    col_trim_true <- "orange"
    
    # =====================
    # 4. Plot
    # =====================
    p <- ggplot(df, aes(x = x, y = y)) +
      
      # ---- Trim (arancione sotto)
      geom_point(
        data = subset(df, trim_plot == TRUE),
        aes(x = x, y = y),
        shape = 16,
        size  = 0.6,
        color = col_trim_true
      ) +
      
      # ---- FAPs (verde sopra)
      geom_point(
        data = subset(df, faps_plot == TRUE),
        aes(x = x, y = y),
        shape = 16,
        size  = 0.6,
        color = col_faps_true
      ) +
      
      coord_equal() +
      scale_y_reverse() +
      theme_void() +
      ggtitle(name)
    
    # =====================
    # 5. Salvataggio
    # =====================
    out_dir <- "plots_boolean"
    dir.create(out_dir, showWarnings = FALSE)
    
    ggsave(
      filename = file.path(out_dir, paste0(name, "_boolean.png")),
      plot = p,
      width = 6,
      height = 5,
      dpi = 300
    )
  })
  








