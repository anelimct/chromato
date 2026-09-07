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

classer_sub_type <- function(data, 
                             iso_col = "isoprene", 
                             mono_col = "monoterpenes",
                             bvocs_col = NULL,
                             seuil_iso_low = 10,
                             seuil_iso_medium = 30,
                             seuil_mono_low = 2,
                             seuil_mono_medium = 5.1,
                             threshold_iso_emit = 1,
                             threshold_mono_emit = 0.1) {
  
  classer <- function(iso, mono) {
    emit_iso <- !is.na(iso) && iso > threshold_iso_emit
    emit_mono <- !is.na(mono) && mono > threshold_mono_emit
    
    if (!emit_iso && !emit_mono) return("NE")
    if (emit_iso && emit_mono) return("both")
    
    if (emit_iso) {
      niveau_iso <- if (iso < seuil_iso_low) "low"
      else if (iso <= seuil_iso_medium) "medium"
      else "high"
      return(paste0("iso_", niveau_iso))
    }
    
    if (emit_mono) {
      niveau_mono <- if (mono < seuil_mono_low) "low"
      else if (mono <= seuil_mono_medium) "medium"
      else "high"
      return(paste0("mono_", niveau_mono))
    }
  }
  
  data$sub_type <- mapply(classer, 
                          data[[iso_col]], 
                          data[[mono_col]], 
                          SIMPLIFY = TRUE)
  
  # Colonne "type" agrégée : iso / mono / both / NE
  data$type <- dplyr::case_when(
    data$sub_type == "NE"                        ~ "NE",
    data$sub_type == "both"                       ~ "both",
    grepl("^iso_", data$sub_type)                  ~ "iso",
    grepl("^mono_", data$sub_type)                 ~ "mono",
    TRUE                                            ~ NA_character_
  )
  
  # Correction : NA (et non "NE"/"iso"/"mono"/"both") quand il n'y a pas de données BVOC
  if (!is.null(bvocs_col) && bvocs_col %in% names(data)) {
    data$sub_type[data[[bvocs_col]] == 0] <- NA
    data$type[data[[bvocs_col]] == 0] <- NA
  }
  
  return(data)
}



calculate_emission_stats <- function(grid_to_plot, working_file, all_data_mean_EF_taxon) {
  
  all_data_mean_EF_taxon <- all_data_mean_EF_taxon |> mutate(Sum = isoprene + monoterpenes)
  
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
      
      # Statistiques pour Sum
      mean_Sum = mean(Sum, na.rm = TRUE),
      median_Sum = median(Sum, na.rm = TRUE),
      sd_Sum = sd(Sum, na.rm = TRUE),
      
      
      #kurtosis 
      
      kurtosis_sum= moments::kurtosis(Sum, na.rm = TRUE), 
      kurtosis_iso= moments::kurtosis(isoprene, na.rm = TRUE),
      kurtosis_mono= moments::kurtosis(monoterpenes, na.rm = TRUE),
      
      # Nombre total d'espèces dans la grille
      total_species_in_grid = dplyr::n(),
      
      # Proportions par type d'émetteur (calculées sur les espèces avec type non-NA)
      n_type_non_na = sum(!is.na(type)),
      prop_iso = sum(type == "iso", na.rm = TRUE) / n_type_non_na,
      prop_mono = sum(type == "mono", na.rm = TRUE) / n_type_non_na,
      prop_both = sum(type == "both", na.rm = TRUE) / n_type_non_na,
      prop_NE = sum(type == "NE", na.rm = TRUE) / n_type_non_na,
      
      .groups = "drop"
    )
  
  return(emission_stats)
}


# map_emission_stats <- function(emission_stats_df, WOODIV_grid, WOODIV_shape, output_dir = "figures/emission_maps") {
#   
#   # S'assurer que le répertoire de sortie existe
#   if (!dir.exists(output_dir)) {
#     dir.create(output_dir, recursive = TRUE)
#   }
#   
#   # Convertir emission_stats_df en tibble si ce n'est pas déjà le cas
#   emission_stats_df <- as_tibble(emission_stats_df)
#   
#   # Joindre les statistiques avec les géométries des grilles
#   emission_stats_sf <- WOODIV_grid |>
#     dplyr::filter(idgrid %in% emission_stats_df$idgrid) |>
#     dplyr::left_join(emission_stats_df, by = "idgrid") |>
#     sf::st_as_sf()
#   
#   # 1. Carte pour mean_isoprene
#   map_mean_isoprene <- ggplot() +
#     geom_sf(data = WOODIV_shape, fill = "lightgrey", color = NA) +
#     geom_sf(data = emission_stats_sf, aes(fill = mean_isoprene), color = NA) +
#     scale_fill_gradient(low = "white", high = "#03a219", name = "Mean Isoprene (µg g⁻¹ h⁻¹)") +
#     theme_minimal() +
#     labs(title = "Mean Isoprene Emission by Grid",
#          fill = "Mean Isoprene (µg g⁻¹ h⁻¹)") +
#     theme(
#       panel.grid.major = element_blank(),
#       panel.grid.minor = element_blank(),
#       legend.position = "bottom"
#     )
#   
#   # 2. Carte pour median_isoprene
#   map_median_isoprene <- ggplot() +
#     geom_sf(data = WOODIV_shape, fill = "lightgrey", color = NA) +
#     geom_sf(data = emission_stats_sf, aes(fill = median_isoprene), color = NA) +
#     scale_fill_gradient(low = "white", high = "#03a219", name = "Median Isoprene (µg g⁻¹ h⁻¹)") +
#     theme_minimal() +
#     labs(title = "Median Isoprene Emission by Grid",
#          fill = "Median Isoprene (µg g⁻¹ h⁻¹)") +
#     theme(
#       panel.grid.major = element_blank(),
#       panel.grid.minor = element_blank(),
#       legend.position = "bottom"
#     )
#   
#   
#   # 3. Carte pour sd_isoprene
#   map_sd_isoprene <- ggplot() +
#     geom_sf(data = WOODIV_shape, fill = "lightgrey", color = NA) +
#     geom_sf(data = emission_stats_sf, aes(fill = sd_isoprene), color = NA) +
#     scale_fill_gradient(low = "white", high = "#03a219", name = "SD Isoprene (µg g⁻¹ h⁻¹)") +
#     theme_minimal() +
#     labs(title = "SD Isoprene Emission by Grid",
#          fill = "SD Isoprene") +
#     theme(
#       panel.grid.major = element_blank(),
#       panel.grid.minor = element_blank(),
#       legend.position = "bottom"
#     )
#   
#   
#   
#   # 4. Carte pour mean_monoterpenes
#   map_mean_monoterpenes <- ggplot() +
#     geom_sf(data = WOODIV_shape, fill = "lightgrey", color = NA) +
#     geom_sf(data = emission_stats_sf, aes(fill = mean_monoterpenes), color = NA) +
#     scale_fill_gradient(low = "white", high = "#0985e2", name = "Mean Monoterpenes (µg g⁻¹ h⁻¹)") +
#     theme_minimal() +
#     labs(title = "Mean Monoterpenes Emission by Grid",
#          fill = "Mean Monoterpenes (µg g⁻¹ h⁻¹)") +
#     theme(
#       panel.grid.major = element_blank(),
#       panel.grid.minor = element_blank(),
#       legend.position = "bottom"
#     )
#   
#   # 5. Carte pour median_monoterpenes
#   map_median_monoterpenes <- ggplot() +
#     geom_sf(data = WOODIV_shape, fill = "lightgrey", color = NA) +
#     geom_sf(data = emission_stats_sf, aes(fill = median_monoterpenes), color = NA) +
#     scale_fill_gradient(low = "white", high = "#0985e2", name = "Median Monoterpenes (µg g⁻¹ h⁻¹)") +
#     theme_minimal() +
#     labs(title = "Median Monoterpenes Emission by Grid",
#          fill = "Median Monoterpenes (µg g⁻¹ h⁻¹)") +
#     theme(
#       panel.grid.major = element_blank(),
#       panel.grid.minor = element_blank(),
#       legend.position = "bottom"
#     )
#   
#   # 6. Carte pour sd_monoterpenes
#   map_sd_monoterpenes <- ggplot() +
#     geom_sf(data = WOODIV_shape, fill = "lightgrey", color = NA) +
#     geom_sf(data = emission_stats_sf, aes(fill = sd_monoterpenes), color = NA) +
#     scale_fill_gradient(low = "white", high = "#0985e2", name = "SD Monoterpenes (µg g⁻¹ h⁻¹)") +
#     theme_minimal() +
#     labs(title = "SD Monoterpenes Emission by Grid",
#          fill = "SD Monoterpenes") +
#     theme(
#       panel.grid.major = element_blank(),
#       panel.grid.minor = element_blank(),
#       legend.position = "bottom"
#     )
#   
#   
#   # 7. Carte pour mean_Sum
#   map_mean_Sum <- ggplot() +
#     geom_sf(data = WOODIV_shape, fill = "lightgrey", color = NA) +
#     geom_sf(data = emission_stats_sf, aes(fill = mean_Sum), color = NA) +
#     scale_fill_gradient(low = "white", high = "#f16700", name = "Mean Sum isoprenoids (µg g⁻¹ h⁻¹)") +
#     theme_minimal() +
#     labs(title = "Mean Sum isoprenoids Emission by Grid",
#          fill = "Sum isoprenoids") +
#     theme(
#       panel.grid.major = element_blank(),
#       panel.grid.minor = element_blank(),
#       legend.position = "bottom"
#     )
#   
#   # 8. Carte pour prop_iso
#   map_prop_iso <- ggplot() +
#     geom_sf(data = WOODIV_shape, fill = "lightgrey", color = NA) +
#     geom_sf(data = emission_stats_sf, aes(fill = prop_iso), color = NA) +
#     scale_fill_gradient(low = "white", high = "#03a219", name = "Proportion", limits = c(0, 1), labels = scales::percent) +
#     theme_minimal() +
#     labs(title = "Proportion of Isoprene-only Emitters by Grid",
#          fill = "Proportion") +
#     theme(
#       panel.grid.major = element_blank(),
#       panel.grid.minor = element_blank(),
#       legend.position = "bottom"
#     )
#   
#   # 9. Carte pour prop_mono
#   map_prop_mono <- ggplot() +
#     geom_sf(data = WOODIV_shape, fill = "lightgrey", color = NA) +
#     geom_sf(data = emission_stats_sf, aes(fill = prop_mono), color = NA) +
#     scale_fill_gradient(low = "white", high = "#0985e2", name = "Proportion", limits = c(0, 1), labels = scales::percent) +
#     theme_minimal() +
#     labs(title = "Proportion of Monoterpenes-only Emitters by Grid",
#          fill = "Proportion") +
#     theme(
#       panel.grid.major = element_blank(),
#       panel.grid.minor = element_blank(),
#       legend.position = "bottom"
#     )
#   
#   # 10. Carte pour prop_both
#   map_prop_both <- ggplot() +
#     geom_sf(data = WOODIV_shape, fill = "lightgrey", color = NA) +
#     geom_sf(data = emission_stats_sf, aes(fill = prop_both), color = NA) +
#     scale_fill_gradient(low = "white", high = "#c34f70", name = "Proportion", limits = c(0, 1), labels = scales::percent) +
#     theme_minimal() +
#     labs(title = "Proportion of Both-Emitters by Grid",
#          fill = "Proportion") +
#     theme(
#       panel.grid.major = element_blank(),
#       panel.grid.minor = element_blank(),
#       legend.position = "bottom"
#     )
#   
#   # 11. Carte pour prop_NE
#   map_prop_NE <- ggplot() +
#     geom_sf(data = WOODIV_shape, fill = "lightgrey", color = NA) +
#     geom_sf(data = emission_stats_sf, aes(fill = prop_NE), color = NA) +
#     scale_fill_gradient(low = "white", high = "#fce72e", name = "Proportion", limits = c(0, 1), labels = scales::percent) +
#     theme_minimal() +
#     labs(title = "Proportion of Non-Emitters by Grid",
#          fill = "Proportion") +
#     theme(
#       panel.grid.major = element_blank(),
#       panel.grid.minor = element_blank(),
#       legend.position = "bottom"
#     )
#   
#   
#   # Sauvegarder les cartes
#   ggsave(
#     filename = "map_mean_isoprene.png",
#     plot = map_mean_isoprene,
#     path = output_dir,
#     width = 16, height = 16, units = "cm", bg = "white"
#   )
#   
#   ggsave(
#     filename = "map_median_isoprene.png",
#     plot = map_median_isoprene,
#     path = output_dir,
#     width = 16, height = 16, units = "cm", bg = "white"
#   )
#   
#   ggsave(
#     filename = "map_mean_monoterpenes.png",
#     plot = map_mean_monoterpenes,
#     path = output_dir,
#     width = 16, height = 16, units = "cm", bg = "white"
#   )
#   
#   ggsave(
#     filename = "map_median_monoterpenes.png",
#     plot = map_median_monoterpenes,
#     path = output_dir,
#     width = 16, height = 16, units = "cm", bg = "white"
#   )
#   
#   ggsave(
#     filename = "map_sd_monoterpenes.png",
#     plot = map_sd_monoterpenes,
#     path = output_dir,
#     width = 16, height = 16, units = "cm", bg = "white"
#   )
#   
#   ggsave(
#     filename = "map_sd_isoprene.png",
#     plot = map_sd_isoprene,
#     path = output_dir,
#     width = 16, height = 16, units = "cm", bg = "white"
#   )
#   
#   ggsave(
#     filename = "map_mean_Sum.png",
#     plot = map_mean_Sum,
#     path = output_dir,
#     width = 16, height = 16, units = "cm", bg = "white"
#   )
#   
#   ggsave(
#     filename = "map_prop_iso.png",
#     plot = map_prop_iso,
#     path = output_dir,
#     width = 16, height = 16, units = "cm", bg = "white"
#   )
#   
#   ggsave(
#     filename = "map_prop_mono.png",
#     plot = map_prop_mono,
#     path = output_dir,
#     width = 16, height = 16, units = "cm", bg = "white"
#   )
#   
#   ggsave(
#     filename = "map_prop_both.png",
#     plot = map_prop_both,
#     path = output_dir,
#     width = 16, height = 16, units = "cm", bg = "white"
#   )
#   
#   ggsave(
#     filename = "map_prop_NE.png",
#     plot = map_prop_NE,
#     path = output_dir,
#     width = 16, height = 16, units = "cm", bg = "white"
#   )
#   
#   # Retourner les cartes dans une liste pour utilisation interactive
#   return(list(
#     map_mean_isoprene = map_mean_isoprene,
#     map_median_isoprene = map_median_isoprene,
#     map_sd_isoprene = map_sd_isoprene,
#     map_mean_monoterpenes = map_mean_monoterpenes,
#     map_median_monoterpenes = map_median_monoterpenes,
#     map_sd_monoterpenes = map_sd_monoterpenes,
#     map_mean_Sum = map_mean_Sum,
#     map_prop_iso = map_prop_iso,
#     map_prop_mono = map_prop_mono,
#     map_prop_both = map_prop_both,
#     map_prop_NE = map_prop_NE
#   ))
# }


map_emission_stats <- function(emission_stats_df, WOODIV_grid, WOODIV_shape, output_dir = "figures/emission_maps") {
  
  if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)
  
  emission_stats_df <- as_tibble(emission_stats_df)
  
  emission_stats_sf <- WOODIV_grid |>
    dplyr::filter(idgrid %in% emission_stats_df$idgrid) |>
    dplyr::left_join(emission_stats_df, by = "idgrid") |>
    sf::st_as_sf()
  
  # Définition des métriques : nom de variable, couleur, label, si c'est une proportion (0-1)
  metrics_info <- tibble::tribble(
    ~var,                  ~color,     ~label,                              ~is_prop,
    "mean_isoprene",       "#03a219",  "Mean Isoprene (µg g⁻¹ h⁻¹)",        FALSE,
    "median_isoprene",     "#03a219",  "Median Isoprene (µg g⁻¹ h⁻¹)",      FALSE,
    "sd_isoprene",         "#03a219",  "SD Isoprene (µg g⁻¹ h⁻¹)",          FALSE,
    "mean_monoterpenes",   "#0985e2",  "Mean Monoterpenes (µg g⁻¹ h⁻¹)",    FALSE,
    "median_monoterpenes", "#0985e2",  "Median Monoterpenes (µg g⁻¹ h⁻¹)",  FALSE,
    "sd_monoterpenes",     "#0985e2",  "SD Monoterpenes (µg g⁻¹ h⁻¹)",      FALSE,
    "mean_Sum",            "#f16700",  "Mean Sum isoprenoids (µg g⁻¹ h⁻¹)", FALSE,
    "prop_iso",            "#03a219",  "Proportion Isoprene-only",          TRUE,
    "prop_mono",           "#0985e2",  "Proportion Monoterpenes-only",      TRUE,
    "prop_both",           "#c34f70",  "Proportion Both-Emitters",          TRUE,
    "prop_NE",             "#fce72e",  "Proportion Non-Emitters",           TRUE
  )
  
  make_map_with_hist <- function(var, color, label, is_prop) {
    
    map <- ggplot() +
      geom_sf(data = WOODIV_shape, fill = "lightgrey", color = NA) +
      geom_sf(data = emission_stats_sf, aes(fill = .data[[var]]), color = NA) +
      theme_minimal() +
      labs(title = label, fill = label) +
      theme(
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.position = "bottom"
      )
    
    map <- if (is_prop) {
      map + scale_fill_gradient(low = "white", high = color, name = "Proportion",
                                limits = c(0, 1), labels = scales::percent)
    } else {
      map + scale_fill_gradient(low = "white", high = color, name = label)
    }
    
    # Histogramme en médaillon
    hist_plot <- ggplot(emission_stats_sf, aes(x = .data[[var]])) +
      geom_histogram(fill = color, color = "white", bins = 20) +
      theme_void() +
      theme(
        plot.background = element_rect(fill = "white", color = "grey40", linewidth = 0.3),
        plot.margin = margin(2, 4, 2, 2)
      )
    
    # Assemblage : histogramme inséré dans le coin haut-droit de la carte, plus resserré
    map + patchwork::inset_element(
      hist_plot,
      left = 0.72, bottom = 0.72, right = 0.99, top = 0.99
    )
  }
  
  plots_list <- purrr::pmap(metrics_info, function(var, color, label, is_prop) {
    make_map_with_hist(var, color, label, is_prop)
  })
  names(plots_list) <- paste0("map_", metrics_info$var)
  
  purrr::iwalk(plots_list, ~ ggsave(
    filename = paste0(.y, ".png"),
    plot = .x,
    path = output_dir,
    width = 16, height = 16, units = "cm", bg = "white"
  ))
  
  return(plots_list)
}