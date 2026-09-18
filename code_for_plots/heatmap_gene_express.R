library(SpatialExperiment)
bin16_list_ann <- readRDS("~/data/bin16_list_ann.rds")
spe <- bin16_list_ann$sham_b1
geni <- c("Foxo1", "Foxo3", "Smad2", "Smad3")
rowSums(counts(spe)[geni,])

lapply(bin16_list_ann, function(spe){
  somme <- rowSums(counts(spe)[geni,])
})


  lista_matrici <- lapply(bin16_list_ann, function(spe) {
    somme <- rowSums(counts(spe)[geni, , drop = FALSE])
    tessuto <- unique(spe$tissue_type)
    df <- as.data.frame(t(somme))
    df$tissue_type <- tessuto
    return(df)
  })


df_totale <- do.call(rbind, lista_matrici)
tabella_finale <- aggregate(. ~ tissue_type, data = df_totale, FUN = mean)
print(tabella_finale)


# Installa se non la possiedi: install.packages("pheatmap")
library(pheatmap)

# 1. Trasforma la tabella in una matrice numerica pura
# (pheatmap richiede che la variabile di raggruppamento sia nei 'rownames')
matrice_heatmap <- as.matrix(tabella_finale[, -1])
rownames(matrice_heatmap) <- tabella_finale$tissue_type

# 2. Crea l'heatmap
pheatmap(
  matrice_heatmap, 
  cluster_rows = TRUE,         # Raggruppa automaticamente i tessuti simili
  cluster_cols = TRUE,         # Raggruppa i geni con pattern simili
  display_numbers = TRUE,      # Inserisce il valore medio numerico dentro le celle
  number_color = "black",      # Colore del testo dei numeri
  color = colorRampPalette(c("blue", "white", "red"))(100), # Gradiente classico blu-bianco-rosso
  main = "Heatmap dell'Espressione Media dei Geni per Tessuto"
)


###################################################################################

lista_matrici <- lapply(bin16_list_ann, function(spe) {
  espr_normalizzata <- logcounts(spe)[geni, , drop = FALSE]
  medie_spe <- rowMeans(espr_normalizzata, na.rm = TRUE)
  tessuto <- unique(spe$tissue_type)
  df <- as.data.frame(t(medie_spe))
  df$tissue_type <- tessuto
  return(df)
})

df_totale <- do.call(rbind, lista_matrici)
tabella_finale_norm <- aggregate(. ~ tissue_type, data = df_totale, FUN = mean)

print(tabella_finale_norm)

library(pheatmap)

matrice_heatmap <- as.matrix(tabella_finale_norm[, -1])
rownames(matrice_heatmap) <- tabella_finale_norm$tissue_type

pheatmap(
  matrice_heatmap, 
  cluster_rows = TRUE,         
  cluster_cols = TRUE,         
  scale = "column",             # <--- IMPORTANTE: calcola lo Z-score per colonna (gene)
  color = colorRampPalette(c("blue", "white", "red"))(100),
  main = "Heatmap su Logcounts (Z-score per Gene)"
)


