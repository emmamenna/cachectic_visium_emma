library(anndataR)
library(SingleCellExperiment)
library(SpatialExperiment)
library(Statial)
library(dplyr)
library(ggplot2)
library(ggsurvfit)
library(survival)
library(tibble)
library(scater)
library(ggspavis)
library(tools)

build_spatial_experiment <- function(fp) {
  sample_name <- basename(fp)
  annd <- anndataR::read_h5ad(fp)
  sce <- SingleCellExperiment::SingleCellExperiment(
    assays = sce$X,
    colData = sce$obs,
    rowData = sce$var,
    ReducedDims = sce$obsm
    )
  return(sce)
}






