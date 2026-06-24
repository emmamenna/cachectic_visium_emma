library(readxl)

nonhd <- read.delim("~/degs_murf_nonhd.tsv")
head(nonhd)
geni_nonhd_up <- rownames(nonhd[nonhd$logFC>0 & nonhd$FDR<0.05,])
geni_nonhd_down <- rownames(nonhd[nonhd$logFC<0 & nonhd$FDR<0.05,])

hd_up <- read_xlsx("~/myo_MURF_treated_untreated_DEGs.xlsx",sheet=3)
geni_hd_up <- hd_up$...1
hd_down <- read_xlsx("~/myo_MURF_treated_untreated_DEGs.xlsx",sheet=4)
geni_hd_down <- hd_down$...1

up_list <- list(
  "Visium" = geni_nonhd_up,
  "HD"  = geni_hd_up
)

down_list <- list(
  "Visium" = geni_nonhd_down,
  "HD"  = geni_hd_down
)

#VENN PLOTS ----------------------------------------------------------------------------------
venn_up <- ggvenn(
  up_list, 
  fill_color = c("#0073C2FF", "#EFC000FF"),
  stroke_size = 0.5, 
  text_size = 6,
  set_name_size = 7
) + 
  ggplot2::ggtitle("Overlap for upregulated DEGs between Murf green and black fibers") +
  ggplot2::theme(
    plot.title = ggplot2::element_text(hjust = 0.5, size = 10, face = "bold") 
  )
print(venn_up)
ggsave(venn_up,file="venn_up.pdf")

venn_down <- ggvenn(
  down_list, 
  fill_color = c("#0073C2FF", "#EFC000FF"),
  stroke_size = 0.5, 
  text_size = 6,
  set_name_size = 7
) + 
  ggplot2::ggtitle("Overlap for downregulated DEGs between Murf green and black fibers") +
  ggplot2::theme(
    plot.title = ggplot2::element_text(hjust = 0.5, size = 10, face = "bold") 
  )
print(venn_down)
ggsave(venn_down,file="venn_down.pdf")

#COMMON GENES --------------------------------------------------------------------------------
common_up <- intersect(geni_hd_up, geni_nonhd_up)
common_down <- intersect(geni_hd_down, geni_nonhd_down)
write_xlsx(
  list(
    Sheet1 = data.frame(DEGs_UP = common_up),
    Sheet2 = data.frame(DEGs_DOWN = common_down)
  ),
  "common_DEGs_murf.xlsx"
)

#LOGFC CORRELATION ----------------------------------------------------------------------------
#prendo tutti i DEGs o solo quelli significativi?
hd_all <- read_xlsx("~/myo_MURF_treated_untreated_DEGs.xlsx",sheet=1)
nonhd$geni <- rownames(nonhd)
geni <- intersect(nonhd$geni, hd_all$SYMBOL)
geni <- c(common_up,common_down)
df <- merge(
  nonhd[, c("geni", "logFC", "FDR")],
  hd_all[, c("SYMBOL", "logFC", "FDR")],
  by.x = "geni",
  by.y = "SYMBOL",
  suffixes = c("_nonHD", "_HD")
)
cols <- ifelse(df$FDR_nonHD < 0.05 & df$FDR_HD < 0.05, "lightgreen",
               ifelse(df$FDR_nonHD < 0.05, "violet",
                      ifelse(df$FDR_HD < 0.05, "yellow", "grey")))
plot(
  df$logFC_nonHD,
  df$logFC_HD,
  col = cols,
  xlab = "logFC visium non HD",
  ylab = "logFC visium HD"
)
abline(a = 0, b = 1, col = "blue", lwd = 1)
legend(
  "topleft",
  legend = c("significativo in entrambi",
             "solo nonHD",
             "solo HD"),
  col = c("lightgreen", "violet", "yellow"),
  pch = 12,
  bty = "n"
)




