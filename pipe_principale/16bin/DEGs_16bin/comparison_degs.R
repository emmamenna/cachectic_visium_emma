library(readxl)
library(writexl)

#Extract multiome DEGS -------------------------------------------------------------------------------
file_name <- "~/per_git/multiome_DEGlist_definitivi.xlsx"
sheet_names <- excel_sheets(file_name)
sheet_names <- sheet_names[! sheet_names %in% c("Myonuclei_IIb", "Myonuclei_IIx", "Myonuclei_IIx_IIa",
                                    "Myonuclei_IIx_IIb", "Myonuclei_MTJ", "Myonuclei_NMJ",
                                    "Myonuclei_Trim63")]

for (sheet in sheet_names) {
  current_data <- read_xlsx(file_name, sheet = sheet)
  genes_up <- current_data$Gene_Name_ID[!is.na(current_data$logFC) & current_data$logFC > 0]
  genes_down <- current_data$Gene_Name_ID[!is.na(current_data$logFC) & current_data$logFC < 0]

  name_up <- paste0("multiome_", sheet, "_up")
  name_down <- paste0("multiome_", sheet, "_down")
  
  assign(name_up, genes_up)
  assign(name_down, genes_down)
}

#SAME FOR DEGS SPATIAL -------------------------------------------------------------------
file_name <- "~/per_git/degs_by_celltype.xlsx"

all_sheets <- excel_sheets(file_name)
total_sheets <- length(all_sheets)

target_indices <- sort(c(seq(3, total_sheets, by = 4), seq(4, total_sheets, by = 4)))
target_indices <- target_indices[target_indices < 17 | target_indices > 44] #tolgo i myo vari
target_sheets <- all_sheets[target_indices]

for (sheet in target_sheets) {
  current_data <- read_xlsx(file_name, sheet = sheet)
  genes <- current_data[[1]]
  genes <- genes[!is.na(genes)]
  vector_name <- paste0("spatial_", sheet)
  assign(vector_name, genes)
}

#Venn diagrams ----------------------------------------------------------------------------
library(ggvenn)

multiome_myo_overall_up <- multiome_Myonuclei_overall_up
multiome_myo_overall_down <- multiome_Myonuclei_overall_down
rm(multiome_Myonuclei_overall_up)
rm(multiome_Myonuclei_overall_down)

# 1. Find all variables in your environment that start with "multiome_"
multiome_vars <- ls(pattern = "^multiome_") #18

pdf("Venn_Diagrams_Automated.pdf", width = 8, height = 6)
lista_intersezioni <- list()
lista_solo_spatial <- list()

for (m_var in multiome_vars) {
  
  # Swap "multiome" for "spatial" to guess the matching spatial variable name
  s_var <- sub("^multiome_", "spatial_", m_var)

  if (exists(s_var)) {
    
    # Estrai i vettori dei geni effettivi
    m_genes <- get(m_var)
    s_genes <- get(s_var)
    
    current_list <- list(
      "Multiome" = m_genes,
      "Spatial"  = s_genes
    )
    
    core_name <- sub("^multiome_", "", m_var)
    nome_foglio <- substr(core_name, 1, 31)
    
    #saving DEGs in common
    geni_comuni <- intersect(m_genes, s_genes)
    if (length(geni_comuni) > 0) {
      df_comuni <- data.frame(Gene_Name = geni_comuni, stringsAsFactors = FALSE)
    } else {
      # Se non ci sono geni in comune, scrive un avviso nel foglio
      df_comuni <- data.frame(Gene_Name = "Nessun gene in comune", stringsAsFactors = FALSE)
    }
    lista_intersezioni[[nome_foglio]] <- df_comuni
    
    geni_esclusivi_spatial <- setdiff(s_genes, m_genes)
    if (length(geni_esclusivi_spatial) > 0) {
      df_solo_spatial <- data.frame(Gene_Name = geni_esclusivi_spatial, stringsAsFactors = FALSE)
    } else {
      df_solo_spatial <- data.frame(Gene_Name = "Nessun gene esclusivo", stringsAsFactors = FALSE)
    }
    lista_solo_spatial[[nome_foglio]] <- df_solo_spatial
    
    # Generate the Venn diagram
    p <- ggvenn(
      current_list, 
      fill_color = c("#0073C2FF", "#EFC000FF"),
      stroke_size = 0.5, 
      text_size = 6,
      set_name_size = 7
    ) + 
      ggplot2::ggtitle(paste("Overlap for:", core_name)) +
      ggplot2::theme(
        plot.title = ggplot2::element_text(hjust = 0.5, size = 18, face = "bold") 
      )
    print(p)
    
  } else {
    message("Skipped: No matching spatial vector found for ", m_var)
  }
  dev.off()
  write_xlsx(lista_intersezioni, path = "Geni_Comuni_Intersezione_Final.xlsx")
  write_xlsx(lista_solo_spatial, path = "Geni_Solo_Spatial_Final.xlsx")
}








