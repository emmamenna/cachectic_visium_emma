load("~/spe_list_8bin.RData")
subset_list <- readRDS("~/bin16_list_ann")
parquet_list <- readRDS("~/barcodes_mapping.rds")
library(SpatialExperiment)
library(scater)
library(ggplot2)
library(ggspavis)
source("~/per_git/inflamm_function.R")

# --- BLOCCO 1 ---
c26foxO_b1 <- infl_transfer(parquet_list[[1]], spe_list$blocco1_c26foxO, subset_list$c26foxO_b1)
#plotCoords(c26foxO_b1)
subset_list$c26foxO_b1 <- c26foxO_b1

sham_b1 <- infl_transfer(parquet_list[[1]], spe_list$blocco1_sham, subset_list$sham_b1)
#plotCoords(sham_b1)
subset_list$sham_b1 <- sham_b1

c26STAT3_b1 <- infl_transfer(parquet_list[[1]], spe_list$blocco1_c26STAT3, subset_list$c26STAT3_b1)
#plotCoords(c26STAT3_b1)
subset_list$c26STAT3_b1 <- c26STAT3_b1

# --- BLOCCO 2 --- 
# DIFFERENT SOURCE
spe_list_016um_b2 <- readRDS("~/spe_list_016um_b2.rds")

idx <- match(colnames(blocco2_list$c26murf1_b2), colnames(spe_list_016um_b2$blocco2_c26murf1))
blocco2_list$c26murf1_b2$to_discard <- spe_list_016um_b2$blocco2_c26murf1$to_discard[idx]
spe <- blocco2_list$c26murf1_b2
plotCoords(spe, annotate = "to_discard")

idx <- match(colnames(blocco2_list$c26SMAD23_b2), colnames(spe_list_016um_b2$blocco2_c26SMAD23))
blocco2_list$c26SMAD23_b2$to_discard <- spe_list_016um_b2$blocco2_c26SMAD23$to_discard[idx]
spe <- blocco2_list$c26SMAD23_b2
plotCoords(spe, annotate = "to_discard")

# idx <- match(colnames(subset_list_updated$c26_b2), colnames(spe_list_016um_b2$blocco2_c26))
# subset_list_updated$c26_b2$to_discard <- spe_list_016um_b2$blocco2_c26$to_discard[idx]
# spe <- subset_list_updated$c26_b2
# plotCoords(spe, annotate = "to_discard")
# 
# subset_list_updated <- lapply(subset_list_updated, function(spe) {
#   keep <- !is.na(spe$to_discard)
#   spe[,keep]
# })
# subset_list_updated$c26murf1_b2 <- 
#   subset_list_updated$c26murf1_b2[,subset_list_updated$c26murf1_b2$to_discard == "FALSE"]
# subset_list_updated$c26SMAD23_b2 <- 
#   subset_list_updated$c26SMAD23_b2[,subset_list_updated$c26SMAD23_b2$to_discard == "FALSE"]
# subset_list_updated$c26_b2 <- 
#   subset_list_updated$c26_b2[,subset_list_updated$c26_b2$to_discard == "FALSE"]

# --- BLOCCO 3 ---
sham_b3 <- infl_transfer(parquet_list[[3]], spe_list$blocco3_sham, subset_list$sham_b3)
#plotCoords(sham_b3)
subset_list$sham_b3 <- sham_b3

c26murf1_b3 <- infl_transfer(parquet_list[[3]], spe_list$blocco3_c26murf1, subset_list$c26murf1_b3)
plotCoords(c26murf1_b3)
subset_list$c26murf1_b3 <- c26murf1_b3

c26STAT3_b3 <- infl_transfer(parquet_list[[3]], spe_list$blocco3_c26STAT3, subset_list$c26STAT3_b3)
#plotCoords(c26STAT3_b3)
subset_list$c26STAT3_b3 <- c26STAT3_b3

# --- BLOCCO 4 ---
c26_b4 <- infl_transfer(parquet_list[[4]], spe_list$blocco4_c26, subset_list$c26_b4)
#plotCoords(c26_b4)
subset_list$c26_b4 <- c26_b4

c26foxO_b4 <- infl_transfer(parquet_list[[4]], spe_list$blocco4_c26foxO, subset_list$c26foxO_b4)
#plotCoords(c26foxO_b4)
subset_list$c26foxO_b4 <- c26foxO_b4

c26SMAD23_b4 <- infl_transfer(parquet_list[[4]], spe_list$blocco4_c26SMAD23, subset_list$c26SMAD23_b4)
#plotCoords(c26SMAD23_b4)
subset_list$c26SMAD23_b4 <- c26SMAD23_b4

# --- BLOCCO 5 ---
c26SMAD23_b5 <- infl_transfer(parquet_list[[5]], spe_list$blocco5_c26SMAD23, subset_list$c26SMAD23_b5)
#plotCoords(c26SMAD23_b5)
subset_list$c26SMAD23_b5 <- c26SMAD23_b5

c26STAT3_b5 <- infl_transfer(parquet_list[[5]], spe_list$blocco5_c26STAT3, subset_list$c26STAT3_b5)
#plotCoords(c26STAT3_b5)
subset_list$c26STAT3_b5 <- c26STAT3_b5

c26murf1_b5 <- infl_transfer(parquet_list[[5]], spe_list$blocco5_c26murf1, subset_list$c26murf1_b5)
#plotCoords(c26murf1_b5)
subset_list$c26murf1_b5 <- c26murf1_b5

# --- BLOCCO 6 ---
c26_b6 <- infl_transfer(parquet_list[[6]], spe_list$blocco6_c26, subset_list$c26_b6)
#plotCoords(c26_b6)
subset_list$c26_b6 <- c26_b6

sham_b6 <- infl_transfer(parquet_list[[6]], spe_list$blocco6_sham, subset_list$sham_b6)
#plotCoords(sham_b6)
subset_list$sham_b6 <- sham_b6

c26foxO_b6 <- infl_transfer(parquet_list[[6]], spe_list$blocco6_c26foxO, subset_list$c26foxO_b6)
plotCoords(c26foxO_b6)
subset_list$c26foxO_b6 <- c26foxO_b6

# --- BLOCCO 9 ---
c26SMAD23_b9 <- infl_transfer(parquet_list[[7]], spe_list$blocco9_c26SMAD23, subset_list$c26SMAD23_b9)
plotCoords(c26SMAD23_b9)
subset_list$c26SMAD23_b9 <- c26SMAD23_b9
c26foxO_b9 <- infl_transfer(parquet_list[[7]], spe_list$blocco9_c26foxO, subset_list$c26foxO_b9)
plotCoords(c26foxO_b9)
subset_list$c26foxO_b9 <- c26foxO_b9
c26murf1_b9 <- infl_transfer(parquet_list[[7]], spe_list$blocco9_c26murf1, subset_list$c26murf1_b9)
plotCoords(c26murf1_b9)
subset_list$c26murf1_b9 <- c26murf1_b9

#controllo e salvo
subset_list
saveRDS(blocco2_list, file = "~/blocco2_list.rds")

discard <- is.na(blocco2_list$c26SMAD23_b2$to_discard) | blocco2_list$c26SMAD23_b2$to_discard == "TRUE" 
blocco2_list$c26SMAD23_b2 <- blocco2_list$c26SMAD23_b2[,!discard]

discard <- is.na(blocco2_list$c26murf1_b2$to_discard) | blocco2_list$c26murf1_b2$to_discard == "TRUE" 
blocco2_list$c26murf1_b2 <- blocco2_list$c26murf1_b2[,!discard]

# lapply(subset_list, function(spe){
#   print(dim(spe))
# })
# 
# results <- lapply(subset_list, function(spe){
#   total_umi <- sum(colSums(counts(spe)))
#   n_geni <- sum(rowSums(counts(spe)) != 0)
#   
#   return(c(total_umi = total_umi, n_geni = n_geni))
# })
# df <- do.call(rbind, results)
# df <- as.data.frame(df)
# df$sample <- names(subset_list)
# library(openxlsx)
# write.xlsx(df, file = "risultati.xlsx", rowNames = FALSE)







