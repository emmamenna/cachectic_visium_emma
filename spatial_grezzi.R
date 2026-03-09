#Librerie

library(Matrix)
library(jsonlite)
library(SummarizedExperiment)
library(S4Vectors)
library(SpatialExperiment)
library(VisiumIO)
library(SingleCellExperiment)

#Acquisizione blocchi

Visium2HD <- function(bin = "008", blocco_num = "blocco1", directory =  "/mnt/europa/data/sandri/241219_A00626_0902_AHWH77DMXY_3/space_out") {
  # definiamo la directory
  dir <- paste0(directory, "/", blocco_num, "/outs")
  
  spe <- TENxVisiumHD(
    spacerangerOut = dir,
    processing = "raw",
    images = c("hires", "lowres"),
    bin_size = bin,
    jsonFile = "scalefactors_json.json",
    tissuePattern = "tissue_positions.parquet",
    spatialCoordsNames = c("pxl_col_in_fullres", "pxl_row_in_fullres")
  ) |> import()
  
  return(spe)
}

load("~/liste_tissue_cond_016.RData")

processa_blocco <- function(blocco_num, lista_contissue) {
  blocco_nome <- paste0("blocco", blocco_num)
  speHD <- Visium2HD(bin = "016", blocco = blocco_nome)
  
  #Filtro per in/out of tissue (da qupath) e distinguo i 3 campioni
  speHD$tissue_type <- lista_contissue[[blocco_nome]]$exp_condition
  speHD_g <- speHD[, lista_contissue[[blocco_nome]]$in_tissue == TRUE]
  
  assign(paste0("speHD", blocco_num, "_g"), speHD_g, envir = .GlobalEnv)
}

processa_blocco(6, lista_contissue)

save(speHD6_g,file="speHD6_g.RData")

