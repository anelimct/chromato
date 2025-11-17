boxplot_EF <- function (data, tree, field){
  
  data_select <- data |>  dplyr::select(Taxon, spagg, gragg, EF, Compound, PAR_algo, T_algo_K) |> dplyr::mutate(source = "literature") |> rbind(field)
  
  data_iso_mono <- data_select|> dplyr::group_by(Taxon) |> dplyr::filter( length(unique(Compound) ) >= 2) |> dplyr::ungroup()
  ## voir pour quelles espèces on n'a pas les deux coumponds
  data_isoprene <- data_iso_mono |>  dplyr::filter(Compound == "isoprene")
  data_mono <- data_iso_mono |>  dplyr::filter(Compound == "monoterpenes")
  
  gragg_order <- tree$tip.label
  
  data_isoprene$gragg <- factor(data_isoprene$gragg, levels = gragg_order)
  data_mono$gragg <- factor(data_mono$gragg, levels = gragg_order)
  
  p_isoprene <- ggplot(data_isoprene, aes(x = EF, y = gragg, color = gragg)) +
    geom_boxplot() +
    geom_jitter(width = 0.2, alpha = 0.5, aes(color = gragg)) +
    theme_minimal() +
    labs(title = "EF distribution per species for Isoprene",
         x = "EF µg·g⁻¹·h⁻¹",
         y = "Species") +
    theme(legend.position = "none",  axis.text.y = element_text(angle = 20, hjust = 1, size = 6))
  
  ggsave(filename = file.path(here::here("figures"), "isoprene_plot.png"), plot = p_isoprene, width = 10, height = 6)
  
  p_mono <- ggplot(data_mono , aes(x = EF, y = gragg, color = gragg)) +
    geom_boxplot() +
    geom_jitter(width = 0.2, alpha = 0.5, aes(color = gragg)) +
    theme_minimal() +
    labs(title = "EF distribution per species for Monoterpenes",
         x = "EF µg·g⁻¹·h⁻¹",
         y = "Species") +
    theme(legend.position = "none",  axis.text.y = element_text(angle = 20, hjust = 1, size = 6))
  
  
  ggsave(filename = file.path(here::here("figures"), "monoterpenes_plot.png"), plot = p_mono, width = 10, height = 6)
  
  
  
  p_iso_field <- ggplot(data_isoprene |>  
                        dplyr::filter(Taxon %in% unique(field$Taxon)), 
                      aes(x = EF, y = Taxon, color = source)) +
    geom_boxplot(aes(color = NULL), color = "grey50") +  # Boxplot neutre
    geom_jitter(width = 0.2, alpha = 0.7, size = 2) +
    scale_color_manual(values = c("field" = "#b51963", "literature" = "#e3b2c1")) +
    theme_minimal() +
    labs(title = "EF distribution (field data highlighted) for isoprene",
         x = "EF",
         y = "Species") +
    theme(legend.position = "top",
          axis.text.y = element_text(angle = 20, hjust = 1, size = 6))
  
  
  ggsave(filename = file.path(here::here("figures"), "field_highlight_iso.png"), plot = p_iso_field, width = 10, height = 6)
  
  
  p_mono_field <- ggplot(data_mono |>  
                          dplyr::filter(Taxon %in% unique(field$Taxon)), 
                        aes(x = EF, y = Taxon, color = source)) +
    geom_boxplot(aes(color = NULL), color = "grey50") +  # Boxplot neutre
    geom_jitter(width = 0.2, alpha = 0.7, size = 2) +
    scale_color_manual(values = c("field" = "#b51963", "literature" = "#e3b2c1")) +
    theme_minimal() +
    labs(title = "EF distribution (field data highlighted) for monoterpenes",
         x = "EF µg·g⁻¹·h⁻",
         y = "Species") +
    theme(legend.position = "top",
          axis.text.y = element_text(angle = 20, hjust = 1, size = 6))
  
  
  ggsave(filename = file.path(here::here("figures"), "field_highlight_mono.png"), plot = p_mono_field, width = 10, height = 6)
  
  
  
  
  # Préparation des données avec décompte par composé
  species_stats <-
    data_iso_mono |> 
    dplyr::group_by(Taxon) |> 
    dplyr::summarise(
      in_field = any(source == "field"),
      in_literature = any(source == "literature")
    ) |> dplyr::ungroup()
  
  # Histogramme comparé
  p_sp_count <- ggplot(species_stats, aes(x = interaction(in_field, in_literature), 
                            fill = interaction(in_field, in_literature))) +
    geom_bar() +
    scale_x_discrete(labels = c("Nouvelles espèces", "Littérature seulement", "Les deux")) +
    scale_fill_manual(
      values = c("#b51963", "#e3b2c1", "#860a3e"),
      guide = "none"  # Cache la légende car les labels sont sur l'axe x
    ) +
    labs(title = "Répartition des espèces par source de données",
         x = "Catégorie",
         y = "Nombre d'espèces") +
    theme_minimal() +
    geom_text(stat = 'count', aes(label = ..count..), vjust = -0.5) +  # Ajoute les comptes
    theme(axis.text.x = element_text(angle = 45, hjust = 1))  # Incline les labels
  
  
  ggsave(filename = file.path(here::here("figures"), "sp_count.png"), plot = p_sp_count, width = 10, height = 6)
  
  
  writexl::write_xlsx(
    data_iso_mono,
    path = here::here("outputs","data_iso_mono.xlsx" )   # Chemin du fichier de sortie
  )
  
  
  
  return(data_iso_mono)
}



plot_EF_sp <- function (data){
  
  filtered_data <- data |> dplyr::filter(Compound == "monoterpenes") |>  dplyr::group_by(gragg) |>  dplyr::filter("field" %in% source & "literature" %in% source) 
  
  for(gragg_name in unique(filtered_data$gragg)) {
    p <- filtered_data |> 
      dplyr::filter(gragg == gragg_name) |>
      ggplot(aes(x = T_algo_K, y = EF, color = source)) +
      geom_point(size = 3) +
      labs(
        title = paste("Standard emission according to Temperature during sampling", gragg_name),
        x = "Temperature (K)",
        y = "EF µg·g⁻¹·h⁻",
        color = "Source"
      ) +
      scale_color_manual(values = c("field" = "#b51963", "literature" = "#e3b2c1")) 
    
    # Sauvegarder le graphique
    
    ggsave(filename = file.path(here::here("figures", "emission_litt_vs_field"), paste0("graph_", gragg_name, ".png")), plot =  p, width = 8, height = 6)
  
  }
  return(data)
}