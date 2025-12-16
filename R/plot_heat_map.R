create_bubble_heatmap <- function(table, log_relative_prop = FALSE, title = "Carte de chaleur des émissions (échelle logarithmique)") {
  
  # 1. Transformer les données
  table_long <- table %>%
    tidyr::pivot_longer(
      cols = -compound,
      names_to = "species",
      values_to = "emission_rate"
    ) %>%
    dplyr::mutate(
      # Convertir en numérique
      emission_rate_numeric = as.numeric(emission_rate),
      # Remplacer NA par 0
      emission_rate_numeric = ifelse(is.na(emission_rate_numeric), 0, emission_rate_numeric)
    )
  
  # 2. Calculer la proportion relative en pourcentage
  table_long <- table_long %>%
    dplyr::group_by(species) %>%
    dplyr::mutate(
      total_species = sum(emission_rate_numeric, na.rm = TRUE),
      relative_proportion = ifelse(total_species > 0,
                                   (emission_rate_numeric / total_species) * 100,
                                   0)
    ) %>%
    dplyr::ungroup()
  
  # 3. Filtrer pour enlever les zéros (pas de bulle pour 0)
  table_filtered <- table_long %>%
    dplyr::filter(emission_rate_numeric > 0)
  
  # Vérifier s'il y a des données
  if(nrow(table_filtered) == 0) {
    warning("Aucune donnée d'émission positive trouvée")
    return(NULL)
  }
  
  # 4. Appliquer log à la proportion relative si demandé
  if(log_relative_prop) {
    table_filtered <- table_filtered %>%
      dplyr::mutate(
        relative_proportion = log10(relative_proportion + 1)
      )
  }
  
  # 5. Créer la bubble heatmap
  p <- ggplot2::ggplot(
    table_filtered,
    ggplot2::aes(
      x = species,
      y = compound,
      size = log10(emission_rate_numeric + 1),
      color = relative_proportion
    )
  ) +
    ggplot2::geom_point(alpha = 0.8) +
    
    # Échelle de taille basée sur log10
    ggplot2::scale_size_continuous(
      name = "Taux d'émission (log10)",
      range = c(1, 12),
      breaks = c(0, 1, 2, 3),
      labels = c("1", "10", "100", "1000")
    ) +
    
    # Échelle de couleurs viridis inferno
    ggplot2::scale_color_viridis_c(
      name = ifelse(log_relative_prop, 
                    "Proportion relative (log10 %)", 
                    "Proportion relative (%)"),
      option = "inferno",
      direction = -1,
      limits = if(!log_relative_prop) c(0, 100) else NULL,
      breaks = if(!log_relative_prop) c(0, 25, 50, 75, 100) else NULL
    ) +
    
    ggplot2::theme_minimal() +
    ggplot2::theme(
      axis.text.x = ggplot2::element_text(angle = 45, hjust = 1, size = 8),
      axis.text.y = ggplot2::element_text(size = 7),
      panel.grid.major = ggplot2::element_line(color = "gray90"),
      legend.position = "right",
      legend.box = "vertical",
      plot.title = ggplot2::element_text(hjust = 0.5)
    ) +
    ggplot2::labs(
      x = "Espèces",
      y = "Composés volatils",
      title = title
    )
  
  return(p)
}

rename_mean_columns <- function(data, sp_screening) {
  name_mapping <- sp_screening %>%
    dplyr::select(spcode, full_scientific_name) %>%
    dplyr::distinct() %>%
    dplyr::mutate(
      clean_name = gsub("_", " ", full_scientific_name),
      old_col_name = paste0("mean_", tolower(spcode))
    )
  
  rename_vector <- setNames(name_mapping$clean_name, name_mapping$old_col_name)
  
  for(i in seq_along(colnames(data))) {
    col_name <- colnames(data)[i]
    if(col_name %in% names(rename_vector)) {
      colnames(data)[i] <- rename_vector[[col_name]]
    }
  }
  
  return(data)
}
