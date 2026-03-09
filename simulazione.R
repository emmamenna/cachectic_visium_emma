load("~/subset_list.RData")
load("~/spe_list_8bin_ann2.RData")
sham_8bin <- spe_list$sham_b1
sham_16bin <- subset_list$sham_b1
multiome_sce <- readRDS("~/multiome_sce.rds")
library(SpatialExperiment)
library(SingleCellExperiment)
library(spacexr)

#check distribution of counts in 16bin for genes that have null counts in 8bin
# counts_16bin_genes_zero_in_8bin <- counts(sham_16bin)[rowMeans(counts(sham_8bin)) == 0,]
# apply(counts_16bin_genes_zero_in_8bin,1,summary)
#sono tutti zeri 

means16 <- rowMeans(counts(sham_16bin))
means8 <- rowMeans(counts(sham_8bin))
common_genes <- common_genes <- intersect(
  rownames(sham_16bin)[means16 != 0],
  rownames(sham_8bin)[means8 != 0]
)
counts_ratio <- means16[common_genes]/means8[common_genes]
counts_ratio <- ifelse(is.na(counts_ratio),0,counts_ratio)
#plot(density(counts_ratio))
var_ratios <- var(counts_ratio)
sd_ratios <- sqrt(var_ratios) #0.36

library(purrr)

simulate_one <- function(i) {
  repeat {
    sim_ratio <- rnorm(1, mean = 4, sd = sd_ratios)
    if(sim_ratio > 0) break
  }
  sim_counts <- round(counts(sham_16bin)/sim_ratio)
  sim_spe <- SpatialExperiment(
    assays = list(counts = sim_counts),
    spatialCoords = spatialCoords(sham_16bin),
    colData = colData(sham_16bin)
  )
  # Salva i geni che vanno a zero
  genes_zero <- rownames(counts(sim_spe)[rowMeans(counts(sim_spe)) == 0,])
  
  rctd_data <- createRctd(
    spatial_experiment = sim_spe,
    reference_experiment = multiome_sce,
    cell_type_col = "subclusters",
    gene_obs_min = 1
  )
  results <- runRctd(rctd_data, rctd_mode = "multi")
  ws <- assay(results, "weights")
  ws <- data.frame(t(as.matrix(ws)))
  colData(sim_spe)[names(ws)] <- ws[colnames(sim_spe), ]
  ids <- names(ws)[apply(ws, 1, which.max)]
  ids <- gsub("\\.([A-z])", " \\1", ids)
  idx <- match(colnames(sim_spe), rownames(ws))
  sim_spe$cell_type <- factor(ids[idx])
  tab <- table(factor(sim_spe$cell_type))
  my_prop <- prop.table(tab)
  
  list(prop = my_prop, null_genes = genes_zero)
}

# Applica la funzione 20 volte
results_list <- map(1:20, simulate_one)

# Estrai le due liste dai risultati
prop <- map(results_list, "prop")
null_genes <- map(results_list, "null_genes")

#ECCO LA SPIEGAZIONEEEE
lapply(null_genes, function(x) length(x))
lapply(prop, function(x) x[1] < 0.3)
#i geni nulli passano da 22012 a 28102 quando la proporzione
#di cellule endoteliali "scavalla" = secondo picco della distribuzione

#vediamo quali sono
geni_nulli_poche_end <- null_genes[[1]]
geni_nulli_tante_end <- null_genes[[5]]
geni_diff <- geni_nulli_tante_end[!(geni_nulli_tante_end %in% geni_nulli_poche_end)]
print(geni_diff)
marker_genes <- readRDS("~/marker_genes_list_singleR_annot.rds")
markers <- unlist(marker_genes)
geni_diff[geni_diff %in% markers]

# sham_8bin <- sham_8bin[,!is.na(sham_8bin$cell_type)]
# prop.8bin <- mean(sham_8bin$cell_type == "Endothelial") #0.57

library(readxl)
profiles <- read_xlsx("~/profiles_sham_b1.xlsx")
library(dplyr)
profiles_end <- profiles %>%
  arrange(desc(Endothelial))  %>%
  head(50)
geni_diff[geni_diff %in% profiles_end$gene]

profiles_myo2b <- profiles %>%
  arrange(desc(Myonuclei_IIb))  %>%
  head(100)
geni_diff[geni_diff %in% profiles_myo2b$gene]

profiles_faps <- profiles %>%
  arrange(desc(FAPs))  %>%
  head(100)
geni_diff[geni_diff %in% profiles_faps$gene]

profiles_imm <- profiles %>%
  arrange(desc(Immune_Cells))  %>%
  head(100)
geni_diff[geni_diff %in% profiles_imm$gene]

