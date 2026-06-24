subset_list <- readRDS("~/bin16_list_ann.rds")
c26_list <- subset_list[names(subset_list) %in% c("sham_b1","c26_b2","sham_b3","c26_b4","c26_b6",
                                                  "sham_b6")]
prop_df <- lapply(names(c26_list), function(name) {
  spe <- c26_list[[name]]
  tab <- table(spe$new_cell_type)
  prop <- prop.table(tab)
  data.frame(
    sample = name,
    cell_type = names(prop),
    proportion = as.numeric(prop)
  )
}) |> bind_rows()

ggplot(prop_df, aes(x = sample, y = proportion, fill = cell_type)) +
  geom_bar(stat = "identity") +
  scale_fill_manual(values = color_vector) +
  theme_minimal() +
  labs(
    title = "Proporzioni dei tipi cellulari per campione",
    x = "Campione (c26_list)",
    y = "Proporzione",
    fill = "Tipo Cellulare"
  ) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

ggsave("cell_type_proportions.png", width = 10, height = 6, dpi = 300)

prop_df_condition <- lapply(c26_list, function(spe) {
  data.frame(
    tissue_type = unique(spe$tissue_type),
    cell_type = spe$new_cell_type
  )
}) |> bind_rows() |> filter(!is.na(cell_type))

prop_df_condition <- prop_df_condition %>%
  dplyr::count(tissue_type, cell_type) %>%
  dplyr::group_by(tissue_type) %>%
  dplyr::mutate(proportion = n / sum(n)) %>%
  ungroup()

ggplot(prop_df_condition, aes(x = tissue_type, y = proportion, fill = cell_type)) +
  geom_bar(stat = "identity") +
  scale_fill_manual(values = color_vector) +
  theme_minimal() +
  labs(
    title = "Proporzioni dei tipi cellulari per condizione",
    x = "Condizione (tissue_type)",
    y = "Proporzione",
    fill = "Tipo Cellulare"
  ) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

ggsave("cell_type_proportions_by_condition.png", width = 10, height = 6, dpi = 300)

prop_per_sample <- purrr::imap_dfr(c26_list, function(spe, name) {
  cd <- as.data.frame(SummarizedExperiment::colData(spe))
  cd |>
    filter(!is.na(new_cell_type)) |>
    dplyr::count(new_cell_type, name = "count") |>
    mutate(proportion  = count / sum(count),
           sample      = name,
           tissue_type = paste(unique(as.character(cd$tissue_type)), collapse = ", ")) |>
    dplyr::rename(cell_type = new_cell_type)
})
prop_per_sample$tissue_type <- factor(
  prop_per_sample$tissue_type,
  levels = c("sham", "c26")
)
ggplot(prop_per_sample, aes(x = tissue_type, y = proportion, fill = tissue_type)) +
  geom_boxplot(outlier.shape = 21) +
  facet_wrap(~ cell_type, scales = "free_y") +
  scale_fill_manual(values = c("sham" = "#1B9E77", "c26" = "#D95F02")) +
  theme_minimal(base_size = 12) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(title = "Distribuzione proporzioni cell_type per tissue_type",
       x = "Tissue type",
       y = "Proporzione")
ggsave("boxplot_cell_types.png", width = 10, height = 6, dpi = 300)

for (spe_name in names(c26_list)) {
  spe <- subset_list[[spe_name]]
  stat1 <- plotCoords(spe, annotate = "new_cell_type", point_size = 0.6) +
    scale_color_manual(values = color_vector) +
    theme(
      legend.key.width = unit(0.5, "lines"),
      legend.key.height = unit(1, "lines")
    )
  ggsave(paste0(spe_name,"_dec.png"), stat1, width = 10, height = 8, dpi = 300)
}

