load("~/subset_list.RData")
load("~/spe_list_8bin_ann2.RData")
sham_8bin <- spe_list$blocco1_sham
sham_16bin <- subset_list$sham_b1
counts_ratio <- rowMeans(counts(sham_16bin))/rowMeans(counts(sham_8bin))
counts_ratio <- ifelse(is.na(counts_ratio),0,counts_ratio)
plot(density(counts_ratio))

sim_counts <- round(counts(sham_16bin)/4) #sample(counts_ratio)
library(SpatialExperiment)
library(SingleCellExperiment)
sim_spe <- SpatialExperiment(
  assays = list(counts = sim_counts),
  spatialCoords = spatialCoords(sham_16bin),
  colData = colData(sham_16bin),
)

#rctd sui dati simulati
multiome_sce <- readRDS("~/multiome_sce.rds")
library(spacexr)
rctd_data <- createRctd(
  spatial_experiment = sim_spe,
  reference_experiment = multiome_sce,
  cell_type_col = "subclusters",
  gene_obs_min = 1
)
#problema perché i geni sono troppo poco espressi
#effettivamente ci sono 1/4 dei bin rispetto all'8bin vero
#comunque funziona con quel parametro rilassato.
results <- runRctd(rctd_data, rctd_mode = "multi")
#weights extraction
ws <- assay(results, "weights")
ws <- data.frame(t(as.matrix(ws)))
colData(sim_spe)[names(ws)] <- ws[colnames(sim_spe), ] 
#Cluster labelling
ids <- names(ws)[apply(ws, 1, which.max)]
ids <- gsub("\\.([A-z])", " \\1", ids)
idx <- match(colnames(sim_spe), rownames(ws))
sim_spe$cell_type <- factor(ids[idx]) 
print(prop.table(table(sim_spe$cell_type)))

#confronto con l'8bin vero
print(prop.table(table(spe_list$sham_b1$cell_type)))




