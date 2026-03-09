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

#DEG  ------------------------------------------------------------------------------
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

#DOTPLOTS --------------------------------------------------------------------------
load("~/markers_deg.RData")
multiome_sce <- readRDS("~/multiome_sce.rds")
library(scater)

for (name in names(markers_deg)) {
  
  genes <- markers_deg[[name]]
  
  p <- plotDots(
    multiome_sham,
    features = genes,
    group = "subclusters",
    assay.type = "logcounts" #provato anche con i counts non cambia molto
  ) + theme(
    axis.text.x = element_text(angle = 45, hjust = 1)
  )
  
  ggsave(
    filename = paste0("dotplot_", name, ".png"),
    plot = p,
    width = 8,
    height = 6
  )
}

#HEATMAP -----------------------------------------------------------------------------
library(ComplexHeatmap)
library(circlize)
multiome_sham <- multiome_sce[,multiome_sce$condition == "SHAM"]

all_genes <- unlist(markers_deg)
names(markers_deg) <- c("m_endo", "m_faps", "m_imm", "m_musc", "m_nersys", "m_peri", "m_sm",
                        "m_teno", "m_m2b", "m_m2x", "m_m2x2a", "m_m2x2b", "m_mmtj", "m_mnmj",
                        "m_trim")
gene_group <- rep(names(markers_deg), lengths(markers_deg))

mat <- logcounts(multiome_sham)[all_genes, ]

# media per subcluster (aggregazione cellule)
group <- multiome_sham$subclusters
mat_avg <- sapply(unique(group), function(g) {
  rowMeans(mat[, group == g, drop = FALSE])
})

# ordina geni per gruppo
ord <- order(gene_group)
mat_avg <- mat_avg[ord, ]
gene_group <- gene_group[ord]

mat_t <- t(mat_avg)

# annotazione ora sulle colonne
col_anno <- HeatmapAnnotation(
  group = gene_group,
  show_annotation_name = FALSE
)

ht <- Heatmap(
  mat_t,
  name = "expr",
  
  top_annotation = col_anno,
  
  cluster_rows = TRUE,
  cluster_columns = FALSE,
  
  show_column_names = FALSE,
  
  column_split = gene_group,
  
  #dimensioni testo
  row_names_gp = grid::gpar(fontsize = 8),
  column_title_gp = grid::gpar(fontsize = 8),
  row_title_gp = grid::gpar(fontsize = 9)
)

png("heatmap_transposed.png",
    width = 3000,
    height = 1500,
    res = 300)

draw(ht)

dev.off()

#HEATMAP CON MARKERS PAPER ------------------------------------------------------------
library(ComplexHeatmap)
library(circlize)
load("~/markers_paper.RData")
multiome_sce <- readRDS("~/multiome_sce.rds")

all_genes <- unlist(markers_paper)
names(markers_paper) <- c("m_endo", "m_faps", "m_imm", "m_musc", "m_peri", "m_sm",
                        "m_teno", "m_m2b", "m_m2x", "m_mmtj", "m_mnmj")

out_genes <- all_genes[!(all_genes %in% rownames(multiome_sce))]
rename_map <- c(
  "Gatsl2" = "Castor2",
  "Mical22" = "Mical2",
  "Fam19a4" = "Tafa4",
  "Map3k/cl" = "Map3k"
)
markers_paper <- lapply(markers_paper, function(x) {
  x[x %in% names(rename_map)] <- rename_map[x[x %in% names(rename_map)]]
  unique(x)
})
all_genes <- all_genes[all_genes %in% rownames(multiome_sce)]
gene_group <- unlist(lapply(names(markers_paper), function(nm) {
  genes <- markers_paper[[nm]]
  genes <- genes[genes %in% all_genes]  # solo quelli presenti
  rep(nm, length(genes))
}))

# estrai matrice espressione
mat <- logcounts(multiome_sce)[all_genes, ]
# media per subcluster (aggregazione cellule)
group <- multiome_sce$subclusters
mat_avg <- sapply(unique(group), function(g) {
  rowMeans(mat[, group == g, drop = FALSE])
})

# ordina geni per gruppo
ord <- order(gene_group)
mat_avg <- mat_avg[ord, ]
gene_group <- gene_group[ord]

mat_t <- t(mat_avg)

# annotazione ora sulle colonne
col_anno <- HeatmapAnnotation(
  group = gene_group,
  show_annotation_name = FALSE
)

ht <- Heatmap(
  mat_t,
  name = "expr",
  
  top_annotation = col_anno,
  
  cluster_rows = TRUE,
  cluster_columns = FALSE,
  
  show_column_names = FALSE,
  
  column_split = gene_group,
  
  # ↓ dimensioni testo
  row_names_gp = grid::gpar(fontsize = 8),
  column_title_gp = grid::gpar(fontsize = 8),
  row_title_gp = grid::gpar(fontsize = 9)
)

png("heatmap_paper_markers.png",
    width = 3000,
    height = 1500,
    res = 300)

draw(ht)

dev.off()

##HEATMAP CON MARCATORI DI IDENTITA' ------------------------------------------------------------
library(ComplexHeatmap)
library(circlize)
markers <- readRDS("~/marker_genes_list_singleR_annot.rds")
multiome_sce <- readRDS("~/multiome_sce.rds")

rename_map <- c(
  "Ncam"   = "Ncam1",
  "Cd140a" = "Pdgfra",
  "Sca1"   = "Atxn1"
)
markers <- lapply(markers, function(x) {
  x[x %in% names(rename_map)] <- rename_map[x[x %in% names(rename_map)]]
  unique(x)
})

all_genes <- unlist(markers)
names(markers) <- c("m_sm", "m_2b", "m_2x", "m_2x2a", "m_m2x2b",
                    "m_mtj", "m_nmj", "m_trim", "m_faps", "m_endo",
                    "m_musc", "m_teno", "m_imm", "m_nersys", "m_peri"
)

gene_group <- rep(names(markers), lengths(markers))
mat <- logcounts(multiome_sce)[all_genes, ]

# media per subcluster (aggregazione cellule)
group <- multiome_sce$subclusters
mat_avg <- sapply(unique(group), function(g) {
  rowMeans(mat[, group == g, drop = FALSE])
})

# ordina geni per gruppo
ord <- order(gene_group)
mat_avg <- mat_avg[ord, ]
gene_group <- gene_group[ord]

mat_t <- t(mat_avg)

# annotazione ora sulle colonne
col_anno <- HeatmapAnnotation(
  group = gene_group,
  show_annotation_name = FALSE
)

ht <- Heatmap(
  mat_t,
  name = "expr",
  
  top_annotation = col_anno,
  
  cluster_rows = TRUE,
  cluster_columns = FALSE,
  
  show_column_names = FALSE,
  
  column_split = gene_group,

  row_names_gp = grid::gpar(fontsize = 8),
  column_title_gp = grid::gpar(fontsize = 8),
  row_title_gp = grid::gpar(fontsize = 9)
)

png("heatmap_identity_markers.png",
    width = 3000,
    height = 1500,
    res = 300)

draw(ht)

dev.off()

#PROIETTO NEL TESSUTO I DEG CHE SEMBRANO MEGLIO DAI DOTPLOT ---------------------------------------
load("~/subset_list.RData")
load("~/markers_deg.RData")

library(ggplot2)
library(SpatialExperiment)
library(patchwork)

spe <- subset_list$sham_b1

for (list_name in names(markers_deg)) {
  
  markers <- markers_deg[[list_name]]
  
  # Prendi al massimo i primi 10
  markers <- markers[1:min(10, length(markers))]
  
  # Tieni solo quelli presenti nello spe
  markers <- markers[markers %in% rownames(spe)]
  
  if (length(markers) == 0) next
  
  expr_matrix <- counts(spe)[markers, , drop = FALSE]
  
  plot_list <- list()
  
  for (i in seq_along(markers)) {
    
    marker <- markers[i]
    expr <- as.numeric(expr_matrix[marker, ])
    
    coords <- as.data.frame(spatialCoords(spe))
    colnames(coords) <- c("x_coord", "y_coord")
    coords$expr <- expr
    
    p <- ggplot(coords, aes(x = x_coord, y = y_coord, color = expr)) +
      geom_point(size = 0.3) +
      coord_fixed() +
      scale_y_reverse() +
      ggtitle(marker) +
      theme_void() +
      theme(plot.title = element_text(size = 8)) +
      scale_color_gradientn(colors = rev(hcl.colors(9, "Rocket")))
    
    plot_list[[i]] <- p
  }
  
  # Dividi in due gruppi (max 5 per riga)
  first_half  <- wrap_plots(plot_list[1:min(5, length(plot_list))], ncol = 5)
  
  if (length(plot_list) > 5) {
    second_half <- wrap_plots(plot_list[6:length(plot_list)], ncol = 5)
    combined_plot <- first_half / second_half
  } else {
    combined_plot <- first_half
  }
  
  ggsave(
    filename = paste0("plot_sham_b1_", list_name, "_top10.png"),
    plot = combined_plot,
    width = 20,
    height = 8,
    dpi = 300
  )
}
