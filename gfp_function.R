gfp_transfer <- function(parquet_file, se8bin, se16bin) {
  # parquet_list: a list, first element must contain $square_008um and $square_016um vectors
  # se8bin: a data.frame with colnames corresponding to 8bin names and a column 'in_treatment'
  # se16bin: a data.frame with colnames corresponding to 16bin names
  
  # Assign to each 8bin the corresponding 16bin label
  se8bin$cell_type16bin <- parquet_file$square_016um[
    match(
      colnames(se8bin),
      parquet_file$square_008um
    )
  ]
  
  # For each 16bin, summarize the treatment status
  in_treat_16 <- tapply(
    se8bin$in_treatment,
    se8bin$cell_type16bin,
    function(x) {
      x <- x[!is.na(x)]
      if (length(x) == 0) {
        NA
      } else if (all(x == 1)) {
        1
      } else if (all(x == 0)) {
        0
      } else {
        2
      }
    }
  )
  
  print(table(in_treat_16)) # for debugging/check
  
  # Assign the summarized in_treatment info to 16bin data.frame columns
  treat_vec <- in_treat_16[colnames(se16bin)]
  se16bin$in_treatment <- as.factor(treat_vec)
  # Remove 'mixed' ("2") bins
  se16bin <- se16bin[, se16bin$in_treatment != "2" & !is.na(se16bin$in_treatment)]
  se16bin$in_treatment <- droplevels(se16bin$in_treatment)
  
  return(se16bin)
}













