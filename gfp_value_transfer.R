#"copying" GFP information from 8 to 16 um bins thanks to barcodes
bin16_list_ann <- readRDS("~/bin16_list_ann.rds")
barcodes_mapping <- readRDS("~/barcodes_mapping.rds")
bin8_list_raw <- readRDS("~/bin8_list_raw.rds")
source("~/per_git/gfp_function.R")
library(ggplot2)
library(ggspavis)

#Blocco1 fox0
c26foxO_b1 <- gfp_value_transfer(barcodes_mapping[[1]], bin8_list_raw$blocco1_c26foxO,
                                 bin16_list_ann$c26foxO_b1)
plotCoords(bin8_list_raw$blocco1_c26foxO, annotate = "GFP_value")
plotCoords(c26foxO_b1, annotate = "gfp_value")
bin16_list_ann$c26foxO_b1 <- c26foxO_b1

c26murf1_b3 <- gfp_value_transfer(barcodes_mapping[[3]], bin8_list_raw$blocco3_c26murf1,
                            bin16_list_ann$c26murf1_b3)
plotCoords(bin8_list_raw$blocco3_c26murf1, annotate = "GFP_value")
plotCoords(c26murf1_b3, annotate = "gfp_value")
bin16_list_ann $c26murf1_b3 <- c26murf1_b3

c26foxO_b4 <- gfp_value_transfer(barcodes_mapping[[4]], bin8_list_raw$blocco4_c26foxO,
                           bin16_list_ann$c26foxO_b4)
plotCoords(bin8_list_raw$blocco4_c26foxO, annotate = "GFP_value")
plotCoords(c26foxO_b4, annotate = "gfp_value")
bin16_list_ann$c26foxO_b4 <- c26foxO_b4

c26SMAD23_b4 <- gfp_value_transfer(barcodes_mapping[[4]], bin8_list_raw$blocco4_c26SMAD23,
                             bin16_list_ann$c26SMAD23_b4)
plotCoords(bin8_list_raw$blocco4_c26SMAD23, annotate = "GFP_value")
plotCoords(c26SMAD23_b4, annotate = "gfp_value")
bin16_list_ann$c26SMAD23_b4 <- c26SMAD23_b4

c26SMAD23_b5 <- gfp_value_transfer(barcodes_mapping[[5]], bin8_list_raw$blocco5_c26SMAD23,
                             bin16_list_ann$c26SMAD23_b5)
plotCoords(bin8_list_raw$blocco5_c26SMAD23, annotate = "GFP_value")
plotCoords(c26SMAD23_b5, annotate = "gfp_value")
bin16_list_ann$c26SMAD23_b5 <- c26SMAD23_b5

c26murf1_b5 <- gfp_value_transfer(barcodes_mapping[[5]], bin8_list_raw$blocco5_c26murf1,
                            bin16_list_ann$c26murf1_b5)
plotCoords(bin8_list_raw$blocco5_c26murf1, annotate = "GFP_value")
plotCoords(c26murf1_b5, annotate = "gfp_value")
bin16_list_ann$c26murf1_b5 <- c26murf1_b5

c26foxO_b6 <- gfp_value_transfer(barcodes_mapping[[6]], bin8_list_raw$blocco6_c26foxO,
                           bin16_list_ann$c26foxO_b6)
plotCoords(bin8_list_raw$blocco6_c26foxO, annotate = "GFP_value")
plotCoords(c26foxO_b6, annotate = "gfp_value")
bin16_list_ann$c26foxO_b6 <- c26foxO_b6

c26foxO_b9 <- gfp_value_transfer(barcodes_mapping[[7]], bin8_list_raw$blocco9_c26foxO,
                           bin16_list_ann $c26foxO_b9)
plotCoords(bin8_list_raw$blocco9_c26foxO, annotate = "GFP_value")
plotCoords(c26foxO_b9, annotate = "gfp_value")
bin16_list_ann$c26foxO_b9 <- c26foxO_b9

c26murf1_b9 <- gfp_value_transfer(barcodes_mapping[[7]], bin8_list_raw$blocco9_c26murf1,
                            bin16_list_ann$c26murf1_b9)
plotCoords(bin8_list_raw$blocco9_c26murf1, annotate = "GFP_value")
plotCoords(c26murf1_b9, annotate = "gfp_value")
bin16_list_ann$c26murf1_b9 <- c26murf1_b9

c26SMAD23_b9 <- gfp_value_transfer(barcodes_mapping[[7]], bin8_list_raw$blocco9_c26SMAD23,
                             bin16_list_ann$c26SMAD23_b9)
plotCoords(bin8_list_raw$blocco9_c26SMAD23, annotate = "GFP_value")
plotCoords(c26SMAD23_b9, annotate = "gfp_value")
bin16_list_ann$c26SMAD23_b9 <- c26SMAD23_b9

saveRDS(bin16_list_ann,"~/bin16_list_ann_gfp.rds")

#QUELLE DEL BLOCCO 2 SONO SBAGLIATE (NON COINCIDONO CON LO 0-1)
# bin16_df <- bin16_list_ann$c26SMAD23_b2
# spe_df <- spe_list_016um_b2$blocco2_c26SMAD23
# bin16_df$gfp_value <- spe_df$GFP_value[match(colnames(bin16_df), colnames(spe_df))]
# plotCoords(spe_df, annotate = "GFP_value")
# plotCoords(bin16_df, annotate = "gfp_value")
# 
# plotCoords(bin16_list_ann$c26SMAD23_b2, annotate = "in_treatment")
# 
# plotCoords(spe_list_016um_b2$blocco2_c26murf1, annotate = "in_treatment")
# plotCoords(spe_list_016um_b2$blocco2_c26murf1, annotate = "GFP_value")




