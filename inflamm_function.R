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

