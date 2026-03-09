load("~/spe_list_8bin.RData")
subset_list_updated <- readRDS("~/subset_list_updated.rds")
load("~/parquet_list.RData")
library(SpatialExperiment)
library(scater)
library(ggplot2)
library(ggspavis)

infl_transfer <- function(parquet_file, se8bin, se16bin) {
  
  # Assign to each 8bin the corresponding 16bin label
  se8bin$cell_type16bin <- parquet_file$square_016um[
    match(
      colnames(se8bin),
      parquet_file$square_008um
    )
  ]

  #gather all 8um bins corresponding to the 16um ones
  all_8um_per_16um <- split(parquet_file$square_008um, parquet_file$square_016um)
  all_8um_per_16um <- lapply(all_8um_per_16um, function(x) unique(x))
  
  #check the ones that are present in 8um bins and only keep a 16um bin where all 4 are not inflamed
  filtered_8um_barcodes <- colnames(se8bin)
  keepers_16um <- names(all_8um_per_16um)[
    sapply(all_8um_per_16um, function(x) all(x %in% filtered_8um_barcodes))
  ]
  
  #subset object
  filtered_spe16bin <- se16bin[, colnames(se16bin) %in% keepers_16um]
}

# --- BLOCCO 1 ---
c26foxO_b1 <- infl_transfer(parquet_list[[1]], spe_list$blocco1_c26foxO, subset_list_updated$c26foxO_b1)
plotCoords(c26foxO_b1)
subset_list_updated$c26foxO_b1 <- c26foxO_b1
sham_b1 <- infl_transfer(parquet_list[[1]], spe_list$blocco1_sham, subset_list_updated$sham_b1)
plotCoords(sham_b1)
subset_list_updated$sham_b1 <- sham_b1
c26STAT3_b1 <- infl_transfer(parquet_list[[1]], spe_list$blocco1_c26STAT3, subset_list_updated$c26STAT3_b1)
plotCoords(c26STAT3_b1)
subset_list_updated$c26STAT3_b1 <- c26STAT3_b1

# --- BLOCCO 2 --- PROBLEMI
# c26SMAD23_b2 <- infl_transfer(parquet_list[[2]], spe_list$blocco2_c26SMAD23, subset_list_updated$c26SMAD23_b2)
# plotCoords(c26SMAD23_b2)
# c26murf1_b2 <- infl_transfer(parquet_list[[2]], spe_list$blocco2_c26murf1, subset_list_updated$c26murf1_b2)
# plotCoords(c26murf1_b2)
# c26_b2 <- infl_transfer(parquet_list[[2]], spe_list$blocco2_c26, subset_list_updated$c26_b2)
# plotCoords(c26_b2)

# --- BLOCCO 3 ---
sham_b3 <- infl_transfer(parquet_list[[3]], spe_list$blocco3_sham, subset_list_updated$sham_b3)
plotCoords(sham_b3)
subset_list_updated$sham_b3 <- sham_b3
c26murf1_b3 <- infl_transfer(parquet_list[[3]], spe_list$blocco3_c26murf1, subset_list_updated$c26murf1_b3)
plotCoords(c26murf1_b3)
subset_list_updated$c26murf1_b3 <- c26murf1_b3
c26STAT3_b3 <- infl_transfer(parquet_list[[3]], spe_list$blocco3_c26STAT3, subset_list_updated$c26STAT3_b3)
plotCoords(c26STAT3_b3)
subset_list_updated$c26STAT3_b3 <- c26STAT3_b3

# --- BLOCCO 4 ---
c26_b4 <- infl_transfer(parquet_list[[4]], spe_list$blocco4_c26, subset_list_updated$c26_b4)
plotCoords(c26_b4)
subset_list_updated$c26_b4 <- c26_b4
c26foxO_b4 <- infl_transfer(parquet_list[[4]], spe_list$blocco4_c26foxO, subset_list_updated$c26foxO_b4)
plotCoords(c26foxO_b4)
subset_list_updated$c26foxO_b4 <- c26foxO_b4
c26SMAD23_b4 <- infl_transfer(parquet_list[[4]], spe_list$blocco4_c26SMAD23, subset_list_updated$c26SMAD23_b4)
plotCoords(c26SMAD23_b4)
subset_list_updated$c26SMAD23_b4 <- c26SMAD23_b4

# --- BLOCCO 5 ---
c26SMAD23_b5 <- infl_transfer(parquet_list[[5]], spe_list$blocco5_c26SMAD23, subset_list_updated$c26SMAD23_b5)
plotCoords(c26SMAD23_b5)
subset_list_updated$c26SMAD23_b5 <- c26SMAD23_b5
c26STAT3_b5 <- infl_transfer(parquet_list[[5]], spe_list$blocco5_c26STAT3, subset_list_updated$c26STAT3_b5)
plotCoords(c26STAT3_b5)
subset_list_updated$c26STAT3_b5 <- c26STAT3_b5
c26murf1_b5 <- infl_transfer(parquet_list[[5]], spe_list$blocco5_c26murf1, subset_list_updated$c26murf1_b5)
plotCoords(c26murf1_b5)
subset_list_updated$c26murf1_b5 <- c26murf1_b5

# --- BLOCCO 6 ---
c26_b6 <- infl_transfer(parquet_list[[6]], spe_list$blocco6_c26, subset_list_updated$c26_b6)
plotCoords(c26_b6)
subset_list_updated$c26_b6 <- c26_b6
sham_b6 <- infl_transfer(parquet_list[[6]], spe_list$blocco6_sham, subset_list_updated$sham_b6)
plotCoords(sham_b6)
subset_list_updated$sham_b6 <- sham_b6
c26foxO_b6 <- infl_transfer(parquet_list[[6]], spe_list$blocco6_c26foxO, subset_list_updated$c26foxO_b6)
plotCoords(c26foxO_b6)
subset_list_updated$c26foxO_b6 <- c26foxO_b6

# --- BLOCCO 9 ---
c26SMAD23_b9 <- infl_transfer(parquet_list[[7]], spe_list$blocco9_c26SMAD23, subset_list_updated$c26SMAD23_b9)
plotCoords(c26SMAD23_b9)
subset_list_updated$c26SMAD23_b9 <- c26SMAD23_b9
c26foxO_b9 <- infl_transfer(parquet_list[[7]], spe_list$blocco9_c26foxO, subset_list_updated$c26foxO_b9)
plotCoords(c26foxO_b9)
subset_list_updated$c26foxO_b9 <- c26foxO_b9
c26murf1_b9 <- infl_transfer(parquet_list[[7]], spe_list$blocco9_c26murf1, subset_list_updated$c26murf1_b9)
plotCoords(c26murf1_b9)
subset_list_updated$c26murf1_b9 <- c26murf1_b9

saveRDS(subset_list_updated, file = "subset_list_updated.rds")






