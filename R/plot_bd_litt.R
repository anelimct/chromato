boxplot_EF <- function (data, tree){
  
  data_iso_mono <- data |> dplyr::group_by(Taxon) |> dplyr::filter( length(unique(Compound) ) >= 2) |> dplyr::ungroup()
  
  data_isoprene <- data_iso_mono |>  dplyr::filter(Compound == "isoprene")
  data_mono <- data_iso_mono |>  dplyr::filter(Compound == "monoterpenes")
  
  taxon_order <- tree$tip.label
  
  data_isoprene$Taxon <- factor(data_isoprene$Taxon, levels = taxon_order)
  data_mono$Taxon <- factor(data_mono$Taxon, levels = taxon_order)
  
  p_isoprene <- ggplot(data_isoprene, aes(x = EF, y = Taxon, color = Taxon)) +
    geom_boxplot() +
    geom_jitter(width = 0.2, alpha = 0.5, aes(color = Taxon)) +
    theme_minimal() +
    labs(title = "EF distribution per species for Isoprene",
         x = "EF",
         y = "Species") +
    theme(legend.position = "none")
  
  ggsave(filename = file.path(here::here("figures"), "isoprene_plot.png"), plot = p_isoprene, width = 10, height = 6)
  
  p_mono <- ggplot(data_mono , aes(x = EF, y = Taxon, color = Taxon)) +
    geom_boxplot() +
    geom_jitter(width = 0.2, alpha = 0.5, aes(color = Taxon)) +
    theme_minimal() +
    labs(title = "EF distribution per species for Monoterpenes",
         x = "EF",
         y = "Species") +
    theme(legend.position = "none")
  
  ggsave(filename = file.path(here::here("figures"), "monoterpenes_plot.png"), plot = p_mono, width = 10, height = 6)
  
  return(data)
}
