#To investigate cell type profiles
spe <- spe_list$blocco1_sham
rctd_data <- createRctd(
  spatial_experiment = spe,
  reference_experiment = sce_ref,
  cell_type_col = "subclusters",
  UMI_min = 15
)

RCTD_obj <- createRctdConfig(rctd_data, max_cores = 4)
RCTD_obj <- fitBulk(RCTD_obj)

profiles_env <- new.env()

fitBulk_custom <- function (RCTD, sample_id, profiles_env) 
{
  bulkData <- spacexr:::prepareBulkData(spacexr:::cell_type_info(RCTD)$info[[1]], 
                                        spacexr:::spatialRNA(RCTD), spacexr:::internal_vars(RCTD)$gene_list_bulk)
  message("fitBulk: decomposing bulk")
  decompose_results <- spacexr:::decompose_full(bulkData$X, sum(spacexr:::nUMI(spacexr:::spatialRNA(RCTD))), 
                                                bulkData$b, verbose = FALSE, constrain = FALSE, MIN_CHANGE = spacexr:::config(RCTD)$MIN_CHANGE_BULK, 
                                                n.iter = 100, bulk_mode = TRUE)
  spacexr:::internal_vars(RCTD)$proportions <- decompose_results$weights
  spacexr:::cell_type_info(RCTD)$renorm <- spacexr:::cell_type_info(RCTD)$info
  profiles <- spacexr:::getNormRef(spacexr:::spatialRNA(RCTD), 
                                   spacexr:::cell_type_info(RCTD)$info[[1]], spacexr:::internal_vars(RCTD)$gene_list_bulk, 
                                   decompose_results$weights)
  
  profiles_env[[sample_id]] <- profiles
  
  spacexr:::cell_type_info(RCTD)$renorm[[1]] <- profiles
  RCTD
}

RCTD_obj <- fitBulk_custom(RCTD_obj, "sham_b1", profiles_env)
library(writexl)
profiles <- as.data.frame(profiles_env[["sham_b1"]])
profiles$gene <- rownames(profiles)
profiles <- profiles[, c("gene", setdiff(colnames(profiles), "gene"))]
#write_xlsx(profiles, "profiles_sham_b1.xlsx")

library(readxl)
profiles <- read_xlsx("~/profiles_sham_b1.xlsx")
markers <- read.table("~/geni_identita.txt")
markers <- as.vector(unlist(markers))
prof_markers <- profiles[profiles$gene %in% markers,]
genes <- prof_markers$gene
library(tidyverse)
std_prof <- t(scale(t(as.matrix(prof_markers[,-1])),center=FALSE))
std_prof <- as.tibble(std_prof)
std_prof$gene <- genes

# medie_riga <- rowMeans(as.matrix(profiles[,-1]))
# raw_std <- rowSds(as.matrix(profiles[,-1]))
# prof_std <- sweep(as.matrix(profiles[,-1]),1,medie_riga,"-")
# prof_std <- sweep(prof_std,1,raw_std,"/")

expr_long <- std_prof %>%
  pivot_longer(
    cols = -gene,
    names_to = "celltype",
    values_to = "expression"
  )

ggplot(expr_long, aes(x = celltype, y = gene)) +
  geom_point(aes(size = expression, color = expression)) +
  scale_color_viridis_c() +
  scale_size_area() + #COMANDO UTILE X LA VISUALIZZAZIONE
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1)
  )

#I geni identità passano i filtri DI RCTD?

markers <- read.table("~/geni_identita.txt")
markers <- as.vector(unlist(markers))
load("~/subset_list.RData")
load("~/spe_list_8bin.RData")
sce_ref <- readRDS("~/multiome_sce.rds")

rctd_16bin <- createRctd(
  spatial_experiment = subset_list$sham_b1,
  reference_experiment = sce_ref,
  cell_type_col = "subclusters"
)

rctd16_markers <- markers %in% rownames(rctd_16bin$spatial_experiment)

rctd_8bin <- createRctd(
  spatial_experiment = spe_list$blocco1_sham,
  reference_experiment = sce_ref,
  cell_type_col = "subclusters",
  UMI_min = 15 #CAMBIATO
)
rctd8_markers <- markers %in% rownames(rctd_8bin$spatial_experiment)

markers_df <- data.frame(
  markers = markers,
  rctd_16 = rctd16_markers,
  rctd_8 = rctd8_markers
)
library(writexl)
write_xlsx(markers_df, path = "rctd_markers.xlsx")

?runRctd

#Check markers proportions

markers <- read.table("~/geni_identita.txt")
markers <- as.vector(unlist(markers))
load("~/subset_list.RData")
load("~/spe_list_8bin.RData")
library(SpatialExperiment)

#8 bin
spe_markers <- spe_list$blocco1_sham[markers,]
markers_depth <- rowSums(counts(spe_markers))
depth <- sum(markers_depth)
markers_prop <- markers_depth/depth *100

#16bin
spe_markers16 <- subset_list$sham_b1[markers,]
markers_depth16 <- rowSums(counts(spe_markers16))
depth16 <- sum(markers_depth16)
markers_prop16 <- markers_depth16/depth16 * 100

df_markers <- data.frame(
  markers = markers,
  bin8 = markers_prop,
  bin16 = markers_prop16
)
library(writexl)
write_xlsx(df_markers, path="markers_percentages.xlsx")

############################ on 16um bins results ###########################################
prova <- results$cell_type_list
results$conf_list
results$min_score
colData(results)$num_cell_types <- sapply(
  colData(results)$cell_type_list,
  length
)

weights <- assay(results, "weights")

# 2. For each possible cell type, extract weights for each cell if that type is listed in the cell_type_list
# We'll iterate through all possible cell types (appear in assay rownames)
for (ct in rownames(weights)) {
  # For each cell, if that cell type is present in its cell_type_list, record the corresponding weight, otherwise NA
  colData(results)[[paste0("weight_", ct)]] <-
    sapply(seq_along(colData(results)$cell_type_list), function(i) {
      if (ct %in% colData(results)$cell_type_list[[i]]) {
        weights[ct, i]
      } else {
        NA
      }
    })
}

cell_type_weights <- lapply(
  seq_along(colData(results)$cell_type_list),
  function(i) {
    celltypes <- colData(results)$cell_type_list[[i]]
    setNames(weights[celltypes, i], celltypes)
  }
)
colData(results)$cell_type_weights <- cell_type_weights

# Compute the sum for each cell
weight_sums <- sapply(colData(results)$cell_type_weights, sum)

# Option 1: Logical vector--Does each cell's sum == 1?
all_equal_one <- abs(weight_sums - 1) < 1e-6  # Use tolerance for floating-point errors
table(weight_sums == 0) #lessgo

zero_weight_cells <- which(weight_sums == 0)  # Indices of such cells
conf_list_for_zero_weights <- colData(results)$conf_list[zero_weight_cells]
conf_list_for_zero_weights
table(unlist(conf_list_for_zero_weights))
#ottimo, tutti cell types FALSE -> PESI ZERO

# conf_list = colData(results)$conf_list
matches <- sapply(
  colData(results)$conf_list,
  function(x) {
    length(x) > 1 && is.logical(x) && !is.na(x[1]) && x[1] == FALSE && any(x[-1] == TRUE, na.rm = TRUE)
  }
)

# The indices (or spot/cell names) where the condition is met
matching_indices <- which(matches)
matching_cell_names <- rownames(colData(results))[matching_indices]

# Optional: show, filter, analyze further
results[,matching_cell_names]$cell_type_weights

cd_df <- colData(spe) |> as_tibble()

cd_df |> 
  filter(cell_type == "Endothelial", Endothelial != 0) |> 
  summarise(mean = mean(Endothelial),
            sd = sd(Endothelial),
            min = min(Endothelial),
            max = max(Endothelial),
            mad = mad(Endothelial))
cd_df |> 
  filter(cell_type == "Endothelial") |> 
  summarise(mean = mean(Endothelial),
            sd = sd(Endothelial),
            min = min(Endothelial),
            max = max(Endothelial),
            mad = mad(Endothelial))

library(dplyr)
library(ggplot2)

# proportions of cell types
tab <- table(spe$cell_type)
prop <- prop.table(tab)

prop_df <- data.frame(
  cell_type = names(prop),
  proportion = as.numeric(prop)
)

# single stacked bar
ggplot(prop_df, aes(x = "Sample", y = proportion, fill = cell_type)) +
  geom_bar(stat = "identity", width = 0.6) +
  scale_fill_manual(values = color_vector) +
  theme_minimal() +
  labs(
    title = "Cell type proportions",
    x = "Sample",
    y = "Proportion",
    fill = "Cell Type"
  )
