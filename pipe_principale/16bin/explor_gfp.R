bin16_list_ann <- readRDS("~/bin16_list_ann.rds")
names(bin16_list_ann)
library(openxlsx)

#Fox0, SMAD23
sample_names <- c("c26foxO_b1", "c26SMAD23_b2", "c26foxO_b4", "c26SMAD23_b4", "c26SMAD23_b5",
                  "c26foxO_b6", "c26SMAD23_b9", "c26foxO_b9")
bin16_list_ann$c26foxO_b1$myo <- bin16_list_ann$c26foxO_b1$new_cell_type %in% c("Myonuclei_IIb", "Myonuclei_IIx",
                                               "Myonuclei_IIx_IIa", "Myonuclei_IIx_IIb",
                                               "Myonuclei_MTJ", "Myonuclei_NMJ", "Myonuclei_Trim63")
prop.table(table(bin16_list_ann$c26foxO_b1$in_treatment[bin16_list_ann$c26foxO_b1$myo],
                 bin16_list_ann$c26foxO_b1$new_cell_type[bin16_list_ann$c26foxO_b1$myo]))

myo_cell_types <- c("Myonuclei_IIb", "Myonuclei_IIx",
                    "Myonuclei_IIx_IIa", "Myonuclei_IIx_IIb",
                    "Myonuclei_MTJ", "Myonuclei_NMJ", "Myonuclei_Trim63")

# Create output file
wb <- createWorkbook()

for (name in sample_names) {
  df <- bin16_list_ann[[name]]
  myo <- df$new_cell_type %in% myo_cell_types
  
  # Map in_treatment: 1 -> "treated", 0 -> "untreated"
  trt <- as.character(df$in_treatment[myo])
  trt <- ifelse(trt == "1", "treated",
                ifelse(trt == "0", "untreated", trt))
  
  # Row-wise proportions (rows sum to 1)
  tab_counts <- table(trt, df$new_cell_type[myo])
  prop <- prop.table(tab_counts, margin = 1)
  tab <- as.data.frame.matrix(prop)
  
  addWorksheet(wb, name)
  writeData(wb, name, tab, rowNames = TRUE)
}


# Save workbook
saveWorkbook(wb, "myo_prop_tables.xlsx", overwrite = TRUE)
