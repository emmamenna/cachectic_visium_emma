#Trasferiamo le informazioni sulla GFP dagli 8bin ai 16bin grazie a barcode mapping
subset_list_updated <- readRDS("~/subset_list_updated.rds")
barcodes_mapping <- readRDS("~/barcodes_mapping.rds")
load("~/spe_list_8bin.RData")
source("~/per_git/gfp_function.R")

library(ggplot2)
library(ggspavis)

#Blocco1 fox0
c26foxO_b1 <- gfp_transfer(barcodes_mapping[[1]], spe_list$blocco1_c26foxO, subset_list_updated$c26foxO_b1)
table(c26foxO_b1$in_treatment)
plotCoords(spe_list$blocco1_c26foxO, annotate = "in_treatment")
plotCoords(c26foxO_b1, annotate = "in_treatment")
subset_list_updated$c26foxO_b1 <- c26foxO_b1

#Blocco1 STAT3
subset_list_updated$c26STAT3_b1$in_treatment <- rep(TRUE,dim(subset_list_updated$c26STAT3_b1)[2])

#Blocco2 HA PROBLEMI !!!

#Blocco3 murf1
c26murf1_b3 <- gfp_transfer(barcodes_mapping[[3]], spe_list$blocco3_c26murf1,
                            subset_list$c26murf1_b3)
table(c26murf1_b3$in_treatment)
plotCoords(spe_list$blocco3_c26murf1, annotate = "in_treatment")
plotCoords(c26murf1_b3, annotate = "in_treatment")
subset_list$c26murf1_b3 <- c26murf1_b3

#Blocco3 STAT3
subset_list$c26STAT3_b3$in_treatment <- rep(TRUE,dim(subset_list$c26STAT3_b3)[2])

#Blocco4 fox0
c26foxO_b4 <- gfp_transfer(barcodes_mapping[[4]], spe_list$blocco4_c26foxO,
                           subset_list_updated$c26foxO_b4)
table(c26foxO_b4$in_treatment)
plotCoords(spe_list$blocco4_c26foxO, annotate = "in_treatment")
plotCoords(c26foxO_b4, annotate = "in_treatment")
subset_list_updated$c26foxO_b4 <- c26foxO_b4

#Blocco4 
c26SMAD23_b4 <- gfp_transfer(barcodes_mapping[[4]], spe_list$blocco4_c26SMAD23,
                             subset_list_updated$c26SMAD23_b4)
table(c26SMAD23_b4$in_treatment)
plotCoords(spe_list$blocco4_c26SMAD23, annotate = "in_treatment")
plotCoords(c26SMAD23_b4, annotate = "in_treatment")
subset_list_updated$c26SMAD23_b4 <- c26SMAD23_b4

#Blocco 4 C26
subset_list_updated$c26_b4$in_treatment <- rep(FALSE,dim(subset_list_updated$c26_b4)[2])

#Blocco5 smad
c26SMAD23_b5 <- gfp_transfer(barcodes_mapping[[5]], spe_list$blocco5_c26SMAD23,
                             subset_list_updated$c26SMAD23_b5)
table(c26SMAD23_b5$in_treatment)
plotCoords(spe_list$blocco5_c26SMAD23, annotate = "in_treatment")
plotCoords(c26SMAD23_b5, annotate = "in_treatment")
subset_list_updated$c26SMAD23_b5 <- c26SMAD23_b5

#Blocco5 murf1
c26murf1_b5 <- gfp_transfer(barcodes_mapping[[5]], spe_list$blocco5_c26murf1,
                            subset_list_updated$c26murf1_b5)
table(c26murf1_b5$in_treatment)
plotCoords(spe_list$blocco5_c26murf1, annotate = "in_treatment")
plotCoords(c26murf1_b5, annotate = "in_treatment")
subset_list_updated$c26murf1_b5 <- c26murf1_b5

#Blocco5 stat3
subset_list_updated$c26STAT3_b5$in_treatment <- rep(TRUE,dim(subset_list_updated$c26STAT3_b5)[2])

#Blocco6 fox0
c26foxO_b6 <- gfp_transfer(barcodes_mapping[[6]], spe_list$blocco6_c26foxO,
                           subset_list_updated$c26foxO_b6)
table(c26foxO_b6$in_treatment)
plotCoords(spe_list$blocco6_c26foxO, annotate = "in_treatment")
plotCoords(c26foxO_b6, annotate = "in_treatment")
subset_list_updated$c26foxO_b6 <- c26foxO_b6

#Blocco6 c26
subset_list_updated$c26_b6$in_treatment <- rep(FALSE,dim(subset_list_updated$c26_b6)[2])

#Blocco9 fox0
c26foxO_b9 <- gfp_transfer(barcodes_mapping[[7]], spe_list$blocco9_c26foxO,
                           subset_list_updated$c26foxO_b9)
table(c26foxO_b9$in_treatment)
plotCoords(spe_list$blocco9_c26foxO, annotate = "in_treatment")
plotCoords(c26foxO_b9, annotate = "in_treatment")
subset_list_updated$c26foxO_b9 <- c26foxO_b9

#Blocco9 murf
c26murf1_b9 <- gfp_transfer(barcodes_mapping[[7]], spe_list$blocco9_c26murf1,
                            subset_list_updated$c26murf1_b9)
table(c26murf1_b9$in_treatment)
plotCoords(spe_list$blocco9_c26murf1, annotate = "in_treatment")
plotCoords(c26murf1_b9, annotate = "in_treatment")
subset_list_updated$c26murf1_b9 <- c26murf1_b9

#Blocco9 smad
c26SMAD23_b9 <- gfp_transfer(barcodes_mapping[[7]], spe_list$blocco9_c26SMAD23,
                             subset_list_updated$c26SMAD23_b9)
table(c26SMAD23_b9$in_treatment)
plotCoords(spe_list$blocco9_c26SMAD23, annotate = "in_treatment")
plotCoords(c26SMAD23_b9, annotate = "in_treatment")
subset_list_updated$c26SMAD23_b9 <- c26SMAD23_b9

saveRDS(subset_list_updated, file = "subset_list_updated.rds")

