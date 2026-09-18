library(dplyr)
library(tidyr)
library(stringr)
library(ggplot2)

all_nuclei_list <- readRDS("~/data/nuclei_list_ann.rds")
color_vector <- c(
  "Endothelial" = "#0072B2",       # blu vivo
  "FAPs" = "#009E73",              # verde vivo
  "Immune_Cells" = "#D55E00",      # arancione-rosso
  "MuSC" = "#E69F00",              # arancione
  "Myonuclei_IIb" = "#CC79A7",     # magenta
  "Myonuclei_IIx" = "#56B4E9",     # azzurro chiaro
  "Myonuclei_IIx_IIa" = "#F0E442", # giallo
  "Myonuclei_IIx_IIb" = "#9B59B6", # viola chiaro
  "Myonuclei_MTJ" = "#0099B4",     # turchese
  "Myonuclei_NMJ" = "#DDAA33",     # senape
  "Myonuclei_Trim63" = "#CC3311",  # rosso scuro
  "Nervous_System" = "#44AA99",    # verde acqua
  "Pericyte" = "#AA4499",          # viola/rosa
  "Smooth_Muscular" = "#332288",   # blu-viola scuro
  "Tenocyte" = "#A3E635"           # verde lime
)

nuclei_list <- all_nuclei_list[names(all_nuclei_list) %in% c("blocco1_sham","blocco2_c26","blocco3_sham","blocco4_c26","blocco6_c26","blocco6_sham")]

# 1. Your original data prep
prop_df_condition <- lapply(nuclei_list, function(spe) {
  spe <- spe[, !spe$to_discard]
  data.frame(
    tissue_type = unique(spe$sampletype),
    cell_type = spe$cell_type
  )
}) |> 
  bind_rows() |> 
  dplyr::filter(!is.na(cell_type))

prop_df_condition <- prop_df_condition %>%
  dplyr::count(tissue_type, cell_type) %>%
  dplyr::group_by(tissue_type) %>%
  dplyr::mutate(proportion = n / sum(n)) %>% ungroup()

prop_df_condition$tissue_type <- factor(
  prop_df_condition$tissue_type,
  levels = c("sham", "c26")
)

# 2. Create wide dataframe
prop_wide <- prop_df_condition %>%
  mutate(pct = sprintf("%5.1f%%", proportion * 100)) %>% 
  select(cell_type, tissue_type, pct) %>%
  pivot_wider(names_from = tissue_type, values_from = pct, values_fill = " 0.0%")

# ---------------------------------------------------------
# 3. WIDEN THE PROPORTIONS (Spacing settings)
# ---------------------------------------------------------
cond_names <- levels(prop_df_condition$tissue_type)

# -> TWEAK THIS: Set the spacing between columns
separator <- " | " 

# -> TWEAK THIS: Add extra spaces after the cell type names (e.g., + 10 spaces)
max_len <- max(nchar(as.character(prop_wide$cell_type)))
cell_col_width <- max(max_len + 2, nchar("Cell Type")) 

# Build the labels
prop_wide$label_text <- apply(prop_wide, 1, function(x) {
  padded_cell <- str_pad(x["cell_type"], width = cell_col_width, side = "right")
  paste(c(padded_cell, x[cond_names]), collapse = separator)
})

legend_labels <- setNames(prop_wide$label_text, prop_wide$cell_type)

# Pad the title so it aligns perfectly with the widened cell types
padded_title <- str_pad("Cell Type", width = cell_col_width, side = "right")
legend_title <- paste0(padded_title, separator, paste(cond_names, collapse = separator))

# ---------------------------------------------------------
# 4. SHRINK THE PLOT & RENDER
# ---------------------------------------------------------
ggplot(prop_df_condition, aes(x = tissue_type, y = proportion, fill = cell_type)) +
  # Shrink the bars slightly so they don't look bulky when the plot is narrowed
  geom_bar(stat = "identity", width = 0.85) + 
  scale_fill_manual(values = color_vector, labels = legend_labels) +
  theme_minimal() +
  labs(
    title = "Cell type proportions by condition - nuclei",
    x = "Condition",
    y = "Proportions",
    fill = legend_title
  ) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.text = element_text(family = "mono", size = 10),
    legend.title = element_text(family = "mono", face = "bold", size = 11),
    
    # -> TWEAK THIS: aspect.ratio forces the dimensions of the bar chart.
    # A value > 1 makes it taller and narrower (shrinking the plot width).
    # Try 1.2, 1.5, or 2.0 to see what looks best for your specific data.
    aspect.ratio = 2
  )


# You may need to increase the width in ggsave to accommodate the wider text
ggsave("cell_type_proportions_by_condition.png", width = 12, height = 6, dpi = 300)

