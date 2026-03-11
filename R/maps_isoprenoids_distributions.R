select_grids_completed <- function(completeness_df, min_species, min_completeness) {
  # Filtrer les grilles selon les deux seuils
  filtered_grids <- completeness_df |> 
    dplyr::filter(
      total_species >= min_species,
      bvocs_completeness_all >= min_completeness
    ) |> 
    dplyr::pull(idgrid)
  
  return(filtered_grids)
}



calculate_emission_stats <- function(grid_to_plot, working_file, all_data_mean_EF_taxon) {
  
  # Filtrer les espèces présentes dans les grilles sélectionnées
  species_in_grids <- working_file |> 
    dplyr::filter(idgrid %in% grid_to_plot) |> 
    dplyr::select(idgrid, gragg) |> 
    dplyr::distinct(idgrid, gragg, .keep_all = TRUE)
  
  # Joindre avec les données d'émission
  grid_emissions <- species_in_grids |> 
    dplyr::left_join(all_data_mean_EF_taxon, by = "gragg")
  
  # Calculer les statistiques par grille
  emission_stats <- grid_emissions |> 
    dplyr::group_by(idgrid) |> 
    dplyr::summarise(
      # Statistiques pour l'isoprène
      mean_isoprene = mean(isoprene, na.rm = TRUE),
      median_isoprene = median(isoprene, na.rm = TRUE),
      sd_isoprene = sd(isoprene, na.rm = TRUE),
    
      
      # Statistiques pour les monoterpènes
      mean_monoterpenes = mean(monoterpenes, na.rm = TRUE),
      median_monoterpenes = median(monoterpenes, na.rm = TRUE),
      sd_monoterpenes = sd(monoterpenes, na.rm = TRUE),
      
      # Nombre total d'espèces dans la grille
      total_species_in_grid = dplyr::n(),
      
      .groups = "drop"
    )
  
  return(emission_stats)
}



map_emission_stats <- function(emission_stats_df, WOODIV_grid, WOODIV_shape, output_dir = "figures/emission_maps") {
  
  # S'assurer que le répertoire de sortie existe
  if (!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE)
  }
  
  # Convertir emission_stats_df en tibble si ce n'est pas déjà le cas
  emission_stats_df <- as_tibble(emission_stats_df)
  
  # Joindre les statistiques avec les géométries des grilles
  emission_stats_sf <- WOODIV_grid |>
    dplyr::filter(idgrid %in% emission_stats_df$idgrid) |>
    dplyr::left_join(emission_stats_df, by = "idgrid") |>
    sf::st_as_sf()
  
  # 1. Carte pour mean_isoprene
  map_mean_isoprene <- ggplot() +
    geom_sf(data = WOODIV_shape, fill = "lightgrey", color = NA) +
    geom_sf(data = emission_stats_sf, aes(fill = mean_isoprene), color = NA) +
    scale_fill_gradient(low = "white", high = "#03a219", name = "Mean Isoprene (µg g⁻¹ h⁻¹)") +
    theme_minimal() +
    labs(title = "Mean Isoprene Emission by Grid",
         fill = "Mean Isoprene (µg g⁻¹ h⁻¹)") +
    theme(
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      legend.position = "bottom"
    )
  
  # 2. Carte pour median_isoprene
  map_median_isoprene <- ggplot() +
    geom_sf(data = WOODIV_shape, fill = "lightgrey", color = NA) +
    geom_sf(data = emission_stats_sf, aes(fill = median_isoprene), color = NA) +
    scale_fill_gradient(low = "white", high = "#03a219", name = "Median Isoprene (µg g⁻¹ h⁻¹)") +
    theme_minimal() +
    labs(title = "Median Isoprene Emission by Grid",
         fill = "Median Isoprene (µg g⁻¹ h⁻¹)") +
    theme(
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      legend.position = "bottom"
    )
  
  # 3. Carte pour mean_monoterpenes
  map_mean_monoterpenes <- ggplot() +
    geom_sf(data = WOODIV_shape, fill = "lightgrey", color = NA) +
    geom_sf(data = emission_stats_sf, aes(fill = mean_monoterpenes), color = NA) +
    scale_fill_gradient(low = "white", high = "#0985e2", name = "Mean Monoterpenes (µg g⁻¹ h⁻¹)") +
    theme_minimal() +
    labs(title = "Mean Monoterpenes Emission by Grid",
         fill = "Mean Monoterpenes (µg g⁻¹ h⁻¹)") +
    theme(
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      legend.position = "bottom"
    )
  
  # 4. Carte pour median_monoterpenes
  map_median_monoterpenes <- ggplot() +
    geom_sf(data = WOODIV_shape, fill = "lightgrey", color = NA) +
    geom_sf(data = emission_stats_sf, aes(fill = median_monoterpenes), color = NA) +
    scale_fill_gradient(low = "white", high = "#0985e2", name = "Median Monoterpenes (µg g⁻¹ h⁻¹)") +
    theme_minimal() +
    labs(title = "Median Monoterpenes Emission by Grid",
         fill = "Median Monoterpenes (µg g⁻¹ h⁻¹)") +
    theme(
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      legend.position = "bottom"
    )
  
  # Sauvegarder les cartes
  ggsave(
    filename = "map_mean_isoprene.png",
    plot = map_mean_isoprene,
    path = output_dir,
    width = 16, height = 16, units = "cm", bg = "white"
  )
  
  ggsave(
    filename = "map_median_isoprene.png",
    plot = map_median_isoprene,
    path = output_dir,
    width = 16, height = 16, units = "cm", bg = "white"
  )
  
  ggsave(
    filename = "map_mean_monoterpenes.png",
    plot = map_mean_monoterpenes,
    path = output_dir,
    width = 16, height = 16, units = "cm", bg = "white"
  )
  
  ggsave(
    filename = "map_median_monoterpenes.png",
    plot = map_median_monoterpenes,
    path = output_dir,
    width = 16, height = 16, units = "cm", bg = "white"
  )
  
  # Retourner les cartes dans une liste pour utilisation interactivehttps://www.lelivrescolaire.fr/page/16858280
  return(list(
    map_mean_isoprene = map_mean_isoprene,
    map_median_isoprene = map_median_isoprene,
    map_mean_monoterpenes = map_mean_monoterpenes,
    map_median_monoterpenes = map_median_monoterpenes
  ))
}


