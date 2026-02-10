library(SingleCellExperiment)
library(Signac)
library(dplyr)
library(ggplot2)
library(scDblFinder)
library(scuttle)
library(scran)
library(edgeR)
library(ggrepel)
library(EDASeq)
library(muscat)

multiome_sce <- readRDS("~/multiome_sce.rds")
multiome_sce
multiome_sham <- multiome_sce[,multiome_sce$condition == "SHAM"]

#CON CLUSTER EXPERIMENT ---------------------------------------------------------------
# library(clusterExperiment)
# #trova logcounts o specifica trasformazione sotto
# clust_spe <- ClusterExperiment(object = as.matrix(logcounts(multiome_sham)),
#                                clusters = multiome_sham$subclusters)
# ce<-makeDendrogram(clust_spe) #capire parametri: reduceMethod="var", nDims=500
# 
# plotDendrogram(ce)
# #i DEG sono riferiti ai nodi del dendogramma, quindi devo ritrovare quelli
# #corrispondenti ad ogni cluster
# #sandri dice di fare questa cosa per i myo mentre gli altri di trattarli normalmente
# #non so se davvero serve questa funzione, a questo punto forse posso solo fare dei DEG 
# #più rifiniti sotto..
# plotHeatmap(ce)
# 
# #molto lento 
# bestDendro <- getBestFeatures(ce,contrastType="Dendro",
#                               DEMethod="edgeR",p.value=0.05,number=NROW(clust_spe))
# head(bestDendro)

#DEG classici --------------------------------------------------------------------------------------
pb_mult <- aggregateAcrossCells(multiome_sham, use.assay.type = "counts", 
                                id=DataFrame(label = multiome_sham$subclusters,
                                             sample = multiome_sham$dataset))

pb_mult$group <- as.factor(pb_mult$subclusters)
y <- DGEList(counts(pb_mult), samples = as.data.frame(colData(pb_mult)))
keep <- filterByExpr(y, group = y$samples$group)
y <- y[keep,, keep.lib.sizes=FALSE]
y <- calcNormFactors(y)
design <- model.matrix(~ 0 + group, data = y$samples)
colnames(design) <- levels(y$samples$group)
y <- estimateDisp(y, design)
fit <- glmFit(y, design)

contr.matrix <- makeContrasts(
  
  Endothelial_vs_rest =
    Endothelial - (FAPs + Immune_Cells + MuSC + Myonuclei_IIb +
                     Myonuclei_IIx + Myonuclei_IIx_IIa + Myonuclei_IIx_IIb +
                     Myonuclei_MTJ + Myonuclei_NMJ + Myonuclei_Trim63 +
                     Nervous_System + Pericyte + Smooth_Muscular + Tenocyte) / 14,
  
  FAPs_vs_rest =
    FAPs - (Endothelial + Immune_Cells + MuSC + Myonuclei_IIb +
              Myonuclei_IIx + Myonuclei_IIx_IIa + Myonuclei_IIx_IIb +
              Myonuclei_MTJ + Myonuclei_NMJ + Myonuclei_Trim63 +
              Nervous_System + Pericyte + Smooth_Muscular + Tenocyte) / 14,
  
  Immune_Cells_vs_rest =
    Immune_Cells - (Endothelial + FAPs + MuSC + Myonuclei_IIb +
                      Myonuclei_IIx + Myonuclei_IIx_IIa + Myonuclei_IIx_IIb +
                      Myonuclei_MTJ + Myonuclei_NMJ + Myonuclei_Trim63 +
                      Nervous_System + Pericyte + Smooth_Muscular + Tenocyte) / 14,
  
  MuSC_vs_rest =
    MuSC - (Endothelial + FAPs + Immune_Cells + Myonuclei_IIb +
              Myonuclei_IIx + Myonuclei_IIx_IIa + Myonuclei_IIx_IIb +
              Myonuclei_MTJ + Myonuclei_NMJ + Myonuclei_Trim63 +
              Nervous_System + Pericyte + Smooth_Muscular + Tenocyte) / 14,
  
  Nervous_System_vs_rest =
    Nervous_System - (Endothelial + FAPs + Immune_Cells + MuSC +
                        Myonuclei_IIb + Myonuclei_IIx + Myonuclei_IIx_IIa +
                        Myonuclei_IIx_IIb + Myonuclei_MTJ + Myonuclei_NMJ +
                        Myonuclei_Trim63 + Pericyte + Smooth_Muscular + Tenocyte) / 14,
  
  Pericyte_vs_rest =
    Pericyte - (Endothelial + FAPs + Immune_Cells + MuSC +
                  Myonuclei_IIb + Myonuclei_IIx + Myonuclei_IIx_IIa +
                  Myonuclei_IIx_IIb + Myonuclei_MTJ + Myonuclei_NMJ +
                  Myonuclei_Trim63 + Nervous_System + Smooth_Muscular + Tenocyte) / 14,
  
  Smooth_Muscular_vs_rest =
    Smooth_Muscular - (Endothelial + FAPs + Immune_Cells + MuSC +
                         Myonuclei_IIb + Myonuclei_IIx + Myonuclei_IIx_IIa +
                         Myonuclei_IIx_IIb + Myonuclei_MTJ + Myonuclei_NMJ +
                         Myonuclei_Trim63 + Nervous_System + Pericyte + Tenocyte) / 14,
  
  Tenocyte_vs_rest =
    Tenocyte - (Endothelial + FAPs + Immune_Cells + MuSC +
                  Myonuclei_IIb + Myonuclei_IIx + Myonuclei_IIx_IIa +
                  Myonuclei_IIx_IIb + Myonuclei_MTJ + Myonuclei_NMJ +
                  Myonuclei_Trim63 + Nervous_System + Pericyte + Smooth_Muscular) / 14,
  
  levels = design
)

lrt <- glmLRT(fit, contrast = contr.matrix[, "MuSC_vs_rest"])
res <- topTags(lrt, n = 50)$table
markers_musc <- rownames(res)
markers_musc_expr <- rowMeans(counts(multiome_sham[markers_musc]))
markers_musc <- sort_by(markers_musc,markers_musc_expr)

lrt <- glmLRT(fit, contrast = contr.matrix[, "Endothelial_vs_rest"])
res <- topTags(lrt, n = 50)$table
markers_endo <- rownames(res)
markers_endo_expr <- rowMeans(counts(multiome_sham[markers_endo]))
markers_endo <- sort_by(markers_endo,markers_endo_expr)

lrt <- glmLRT(fit, contrast = contr.matrix[, "Tenocyte_vs_rest"])
res <- topTags(lrt, n = 50)$table
markers_teno <- rownames(res)
markers_teno_expr <- rowMeans(counts(multiome_sham[markers_teno]))
markers_teno <- sort_by(markers_teno,markers_teno_expr)

lrt <- glmLRT(fit, contrast = contr.matrix[, "Smooth_Muscular_vs_rest"])
res <- topTags(lrt, n = 50)$table
markers_sm <- rownames(res)
markers_sm_expr <- rowMeans(counts(multiome_sham[markers_sm]))
markers_sm <- sort_by(markers_sm,markers_sm_expr)

lrt <- glmLRT(fit, contrast = contr.matrix[, "Pericyte_vs_rest"])
res <- topTags(lrt, n = 50)$table
markers_peri <- rownames(res)
markers_peri_expr <- rowMeans(counts(multiome_sham[markers_peri]))
markers_peri <- sort_by(markers_peri,markers_peri_expr)

lrt <- glmLRT(fit, contrast = contr.matrix[, "Nervous_System_vs_rest"])
res <- topTags(lrt, n = 50)$table
markers_nersys <- rownames(res)
markers_nersys_expr <- rowMeans(counts(multiome_sham[markers_nersys]))
markers_nersys <- sort_by(markers_nersys,markers_nersys_expr)

lrt <- glmLRT(fit, contrast = contr.matrix[, "Immune_Cells_vs_rest"])
res <- topTags(lrt, n = 50)$table
markers_imm <- rownames(res)
markers_imm_expr <- rowMeans(counts(multiome_sham[markers_imm]))
markers_imm <- sort_by(markers_imm,markers_imm_expr)

lrt <- glmLRT(fit, contrast = contr.matrix[, "FAPs_vs_rest"])
res <- topTags(lrt, n = 50)$table
markers_faps <- rownames(res)
markers_faps_expr <- rowMeans(counts(multiome_sham[markers_faps]))
markers_faps <- sort_by(markers_faps,markers_faps_expr)

#I mionuclei li confronto coi mionuclei
contr.matrix <- makeContrasts(
  Myonuclei_IIb_vs_rest =
    Myonuclei_IIb - (Myonuclei_IIx + Myonuclei_IIx_IIa + Myonuclei_IIx_IIb +
                     Myonuclei_MTJ + Myonuclei_NMJ + Myonuclei_Trim63) / 6,

  Myonuclei_IIx_vs_rest =
   Myonuclei_IIx - (Myonuclei_IIb + Myonuclei_IIx_IIa + Myonuclei_IIx_IIb +
                     Myonuclei_MTJ + Myonuclei_NMJ + Myonuclei_Trim63) / 6,

  Myonuclei_IIx_IIa_vs_rest =
    Myonuclei_IIx_IIa - (Myonuclei_IIb + Myonuclei_IIx + Myonuclei_IIx_IIb +
                         Myonuclei_MTJ + Myonuclei_NMJ + Myonuclei_Trim63) / 6,

  Myonuclei_IIx_IIb_vs_rest =
    Myonuclei_IIx_IIb - (Myonuclei_IIb + Myonuclei_IIx + Myonuclei_IIx_IIa +
                         Myonuclei_MTJ + Myonuclei_NMJ + Myonuclei_Trim63) / 6,

  Myonuclei_MTJ_vs_rest =
    Myonuclei_MTJ - (Myonuclei_IIb + Myonuclei_IIx + Myonuclei_IIx_IIa +
                     Myonuclei_IIx_IIb + Myonuclei_NMJ + Myonuclei_Trim63) / 6,

  Myonuclei_NMJ_vs_rest =
   Myonuclei_NMJ - (Myonuclei_IIb + Myonuclei_IIx + Myonuclei_IIx_IIa +
                     Myonuclei_IIx_IIb + Myonuclei_MTJ + Myonuclei_Trim63) / 6,
  levels = design
)
lrt <- glmLRT(fit, contrast = contr.matrix[, "Myonuclei_IIb_vs_rest"])
res <- topTags(lrt, n = 50)$table
markers_m2b <- rownames(res)
markers_m2b_expr <- rowMeans(counts(multiome_sham[markers_m2b]))
markers_m2b <- sort_by(markers_m2b,markers_m2b_expr)

lrt <- glmLRT(fit, contrast = contr.matrix[, "Myonuclei_IIx_vs_rest"])
res <- topTags(lrt, n = 50)$table
markers_m2x <- rownames(res)
markers_m2x_expr <- rowMeans(counts(multiome_sham[markers_m2x]))
markers_m2x <- sort_by(markers_m2x,markers_m2x_expr)

lrt <- glmLRT(fit, contrast = contr.matrix[, "Myonuclei_IIx_IIa_vs_rest"])
res <- topTags(lrt, n = 50)$table
markers_m2x2a <- rownames(res)
markers_m2x2a_expr <- rowMeans(counts(multiome_sham[markers_m2x2a]))
markers_m2x2a <- sort_by(markers_m2x2a,markers_m2x2a_expr)

lrt <- glmLRT(fit, contrast = contr.matrix[, "Myonuclei_IIx_IIb_vs_rest"])
res <- topTags(lrt, n = 50)$table
markers_m2x2b <- rownames(res)
markers_m2x2b_expr <- rowMeans(counts(multiome_sham[markers_m2x2b]))
markers_m2x2b <- sort_by(markers_m2x2b,markers_m2x2b_expr)

lrt <- glmLRT(fit, contrast = contr.matrix[, "Myonuclei_MTJ_vs_rest"])
res <- topTags(lrt, n = 50)$table
markers_mmtj <- rownames(res)
markers_mmtj_expr <- rowMeans(counts(multiome_sham[markers_mmtj]))
markers_mmtj <- sort_by(markers_mmtj,markers_mmtj_expr)

lrt <- glmLRT(fit, contrast = contr.matrix[, "Myonuclei_NMJ_vs_rest"])
res <- topTags(lrt, n = 50)$table
markers_mnmj <- rownames(res)
markers_mnmj_expr <- rowMeans(counts(multiome_sham[markers_mnmj]))
markers_mnmj <- sort_by(markers_mnmj,markers_mnmj_expr)

#TRIM dal confronto con gli overall sham -------------------------------------------------
pb_mult <- aggregateAcrossCells(multiome_sce, use.assay.type = "counts", 
                                id=DataFrame(label = multiome_sce$subclusters,
                                             sample = multiome_sce$dataset))
table(multiome_sce$condition)
pb_mult$group <- paste0(pb_mult$subclusters, "_", pb_mult$condition)
mio_ov_sham <- unique(grep(paste0("^", "Myonuclei", ".*", "SHAM"), pb_mult$group, value = TRUE))

pb_sub <- pb_mult[, pb_mult$group %in% c("Myonuclei_Trim63_CAC", mio_ov_sham)]
pb_sub$group[pb_sub$group %in% mio_ov_sham] <- "Myonuclei_SHAM"
pb_sub$group <- as.factor(pb_sub$group)
pb_comb <- aggregateData(pb_sub, assay = "counts", fun = "sum", by = "sample")
colData(pb_comb)$group <- pb_sub$group[match(colnames(pb_comb), pb_sub$sample)]
names(assays(pb_comb)) <- "counts"

y <- DGEList(counts(pb_comb), samples = as.data.frame(colData(pb_comb)))
keep <- filterByExpr(y, group = y$samples$group)
y <- y[keep,, keep.lib.sizes=FALSE]
y <- calcNormFactors(y)
design <- model.matrix(~ group, data=y$samples)
y <- estimateDisp(y, design)
fit <- glmFit(y, design)
lrt <- glmLRT(fit, coef = 2)
res <- topTags(lrt, n = 50)$table
markers_trim <- rownames(res)
multiome_cac <- multiome_sce[,multiome_sce$condition == "CAC"]
markers_trim_expr <- rowMeans(counts(multiome_cac[markers_trim]))
markers_trim <- sort_by(markers_trim,markers_trim_expr)

#SALVO I MARKERS ----------------------------------------------------------------------------
markers_deg <- list(markers_endo,markers_faps,markers_imm,markers_musc,
                    markers_nersys,markers_peri,markers_sm,markers_teno,
                    markers_m2b,markers_m2x,markers_m2x2a,markers_m2x2b,
                    markers_mmtj,markers_mnmj, markers_trim)
names(markers_deg) <- c("markers_endo","markers_faps","markers_imm","markers_musc",
                        "markers_nersys","markers_peri","markers_sm","markers_teno",
                        "markers_m2b","markers_m2x","markers_m2x2a","markers_m2x2b",
                        "markers_mmtj","markers_mnmj","markers_trim")
save(markers_deg,file="markers_deg.RData")

#MARKERS DA PAPER PASCAL MAIRE ----------------------------------------------------------
markers_m2x <- c("Myh1", "Tbc1d1", "2310016D03Rik", "Nrp1", "Gbe1", "Tmem65", "Raver2", "Myom2",
            "Kcnn2", "Nos1", "Ptprk", "Gatsl2", "Padi2", "Zfp385b", "Slc8a1", "Kcnq5", "B3galt1",
            "Aldh1a1", "Pde10a")
markers_m2b <- c("Myh4", "Mybpc2", "Myhas", "Pvalb", "Actn3", "Prkag3", "Mylk4", "Sorbs2", "Gpd2",
            "Tpm1", "Tmem233", "Pde4d", "Pfkfb3", "Mical22", "Atp8a1", "Prr33", "Ampd1",
            "9530026P05Rik", "Enox2", "Adamtsl1")
markers_mmtj <- c("Col22a1", "Maml2", "App", "Clstn2", "Prom1", "Usp6nl", "Map3k/cl",
                  "Med12l", "Gm34907", "Rgcc", "Slc24a2", "Net1", "Ankrd1", "Fras1",
                  "Pdzd2", "Fgf7", "Nav3", "Atp13a5", "Asap2", "Lama2", "Col24a1",
                  "Magi1", "Dennd5b", "Frem2", "Adamts20", "Atp1b4", "Atp10a", "Ltbp1",
                  "Thsd7b", "Ncam1")
markers_mnmj <- c("Ano4", "Phldb2", "Epb41l4a", "Hs6st2", "Colq", "Galnt18", "Gm13387",
                  "B4galnt3", "Col13a1", "Rps6ka1", "Chrne", "Ufsp1", "Hs3st3a1", "Chrna1",
                  "Prkar1a", "Pdzrn4", "Col19a1", "Ctdspl", "Ncam1", "Col4a3", "Pla2g7",
                  "Calcrl", "Fam19a4", "Vav3", "Utrn", "Nav3", "Akap11")
markers_faps <- c("Abca8a", "Plxdc2")
markers_musc <- c("Hs6st3", "Pax7")
markers_teno <- c("Itgbl1", "Thbs4", "Cacna1c")
markers_sm <- c("Cacna1c", "Trpc3")
markers_endo <- c("Pecam1","St6galnac3")
markers_imm <- c("Mctp1", "F13a1", "Pde3b")
markers_faps <- c("Pde3b", "Nnat")
markers_peri <- c("Mmrn1", "Rein")

markers_paper <- list(markers_endo,markers_faps,markers_imm,markers_musc,
                      markers_peri,markers_sm,markers_teno,
                    markers_m2b,markers_m2x,markers_mmtj,markers_mnmj)
names(markers_paper) <- c("markers_endo","markers_faps","markers_imm","markers_musc",
                        "markers_peri","markers_sm","markers_teno",
                        "markers_m2b","markers_m2x", "markers_mmtj","markers_mnmj")

markers_deg_par <- markers_deg[names(markers_deg) %in% names(markers_paper)]
overlap <- Map('%in%', markers_paper, markers_deg_par)
lapply(overlap, table)

##confermo che i myonuclei fanno schifo il resto +- c'è
save(markers_paper, file="markers_paper.RData")


