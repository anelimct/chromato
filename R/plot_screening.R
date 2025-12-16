wide_table_sum_per_sample <- function(compounds_table, bvocs_samples, valid_samples_mono, paradise_reports_mono_ER){
  #samples pour lesquels les mono ont été analysés (.cdf)
  rows_to_plot <- names(fusionner_listes(paradise_reports_mono_ER))[grepl("^[a-zA-Z]{4}_", names(fusionner_listes(paradise_reports_mono_ER)))] 
  #samples qui respectent les condition <43 °C , PAR min (pas.cdf avec des majuscules)
  rows_to_keep <- valid_samples_mono$ID |>  stringr::str_to_lower() |>  stringr::str_c(".cdf")
  
  rows_valid <- valid_samples_mono$ID
  
  df <- compounds_table
  colonnes_samples <- names(df)[grepl("^[a-zA-Z]{4}_", names(df))]
  
  # COMPTAGE DE TOUS LES TYPES DE COMPOUNDS (avant conversion numérique)
  # On compte chaque compound détecté (présent ou avec une valeur) comme 1
  all_counts <- compounds_table |> 
    dplyr::filter(class %in% c("Sesquiterpenes", "Oxygenated-sesquiterpenes", 
                               "Monoterpenes", "Oxygenated-monoterpenes")) |> 
    dplyr::group_by(class) |> 
    dplyr::summarise(
      dplyr::across(
        tidyselect::all_of(colonnes_samples), 
        ~ sum(.x != "" & !is.na(.x)),  # Compter si non-vide et non-NA
        .names = "{.col}"
      )
    ) |> 
    dplyr::ungroup() |> 
    tidyr::pivot_longer(
      cols = -class,
      names_to = "sample",
      values_to = "count"
    ) |> 
    tidyr::pivot_wider(
      names_from = class,
      values_from = count,
      names_prefix = "Count_"
    ) |> 
    dplyr::filter(sample %in% rows_to_keep) |> 
    dplyr::mutate(
      ID = stringr::str_replace(sample, ".cdf", ""),
      ID = stringr::str_to_upper(ID)
    )
  
  ## AGRÉGATION DES VALEURS NUMÉRIQUES pour les autres classes
  df_agg <- compounds_table |> 
    dplyr::group_by(class) |> 
    dplyr::summarise(
      dplyr::across(
        tidyselect::all_of(colonnes_samples), 
        ~ sum(as.numeric(.x), na.rm = TRUE),
        .names = "{.col}"
      )
    ) |> 
    dplyr::ungroup() |> 
    tidyr::pivot_longer(
      cols = -class,
      names_to = "sample",
      values_to = "value"
    ) |> 
    tidyr::pivot_wider(
      names_from = class,
      values_from = value
    ) |> 
    dplyr::filter(sample %in% rows_to_keep) |> 
    dplyr::mutate(
      ID = stringr::str_replace(sample, ".cdf", ""),
      ID = stringr::str_to_upper(ID)
    )
  
  # Fusionner les comptes avec les valeurs agrégées
  df_combined <- df_agg |> 
    dplyr::left_join(all_counts |> dplyr::select(-sample), by = "ID")
  
  # Vérifier que toutes les colonnes de comptage existent, sinon les créer avec 0
  count_cols <- c("Count_Sesquiterpenes", "Count_Oxygenated-sesquiterpenes",
                  "Count_Monoterpenes", "Count_Oxygenated-monoterpenes")
  
  for(col in count_cols) {
    if(!col %in% names(df_combined)) {
      df_combined[[col]] <- 0
    }
  }
  
  # Vérifier que les colonnes nécessaires existent pour les monoterpènes
  # Si elles n'existent pas, les créer avec 0
  if("Monoterpenes" %in% names(df_combined) && "Oxygenated-monoterpenes" %in% names(df_combined)) {
    df_combined <- df_combined |> 
      dplyr::mutate(
        Total_monoterpenes = Monoterpenes + `Oxygenated-monoterpenes`
      )
  } else if("Monoterpenes" %in% names(df_combined)) {
    df_combined <- df_combined |> 
      dplyr::mutate(
        `Oxygenated-monoterpenes` = 0,
        Total_monoterpenes = Monoterpenes
      )
  } else if("Oxygenated-monoterpenes" %in% names(df_combined)) {
    df_combined <- df_combined |> 
      dplyr::mutate(
        Monoterpenes = 0,
        Total_monoterpenes = `Oxygenated-monoterpenes`
      )
  } else {
    df_combined <- df_combined |> 
      dplyr::mutate(
        Monoterpenes = 0,
        `Oxygenated-monoterpenes` = 0,
        Total_monoterpenes = 0
      )
  }
  
  # Calculer le compte total des monoterpènes
  df_combined <- df_combined |> 
    dplyr::mutate(
      Count_Total_monoterpenes = Count_Monoterpenes + `Count_Oxygenated-monoterpenes`
    )
  
  # S'assurer qu'Isoprene existe
  if(!"Isoprene" %in% names(df_combined)) {
    df_combined$Isoprene <- 0
  }
  
  cols_PAR <- c("PAR.début.1", "PAR.début.2", "PAR.milieu.1",     
                "PAR.milieu.2", "PAR.milieu.3", "PAR.milieu.4",  
                "PAR.milieu.5", "PAR.milieu.6", "PAR.fin.1", "PAR.fin.2")
  
  # Format wide final
  bvocs_samples_wide <- bvocs_samples |> 
    dplyr::inner_join(df_combined, by = "ID") |> 
    dplyr::filter(!is.na(Leaves_DM)) |> 
    dplyr::mutate(dplyr::across(dplyr::all_of(cols_PAR), as.numeric)) |> 
    dplyr::mutate(
      PAR_algo = rowMeans(dplyr::across(dplyr::all_of(cols_PAR)), na.rm = TRUE),
      mean_T = purrr::map_dbl(values_T_in, ~ mean(.x, na.rm = TRUE)),
      T_algo_K = mean_T + 273.15,
      Standardized = TRUE
    ) |> 
    dplyr::select(
      ID, 
      Isoprene, 
      Monoterpenes,
      `Oxygenated-monoterpenes`,
      Total_monoterpenes,
      Count_Monoterpenes,
      `Count_Oxygenated-monoterpenes`,
      Count_Total_monoterpenes,
      Count_Sesquiterpenes,
      `Count_Oxygenated-sesquiterpenes`,
      Taxon, 
      PAR_algo, 
      T_algo_K,
      Standardized
    ) |> 
    dplyr::mutate(Taxon = stringr::str_replace_all(Taxon, " ", "_"))
  
  return(bvocs_samples_wide)
}



create_compound_boxplots <- function(data, trait_col, save_plots = FALSE, save_dir = NULL) {
  
  # Palette de 2 couleurs pour les groupes
  color_palette <- c("#bd1b6e", "#ff862f")  # Rose, Orange
  
  # Supprimer les lignes où trait == "intermediate"
  data <- data |> 
    dplyr::filter(tolower(!!sym(trait_col)) != "intermediate")
  
  # Obtenir les groupes uniques
  unique_groups <- unique(data[[trait_col]])
  
  # Garder seulement les 2 premiers groupes
  if(length(unique_groups) > 2) {
    unique_groups <- unique_groups[1:2]
    data <- data |> 
      dplyr::filter(!!sym(trait_col) %in% unique_groups)
  }
  
  group_colors <- setNames(
    color_palette[1:length(unique_groups)],
    unique_groups
  )
  
  # Séparer les colonnes en deux catégories
  emission_cols <- c("Isoprene", "Monoterpenes", "Oxygenated-monoterpenes", "Total_monoterpenes")
  richness_cols <- c("Count_Monoterpenes", "Count_Oxygenated-monoterpenes", "Count_Total_monoterpenes",
                     "Count_Sesquiterpenes", "Count_Oxygenated-sesquiterpenes")
  
  # Filtrer les colonnes disponibles
  available_emission <- emission_cols[emission_cols %in% names(data)]
  available_richness <- richness_cols[richness_cols %in% names(data)]
  
  # Créer les plots d'émissions (avec log)
  emission_plots <- list()
  if(length(available_emission) > 0) {
    emission_plots <- create_emission_boxplots(data, trait_col, available_emission, group_colors)
  }
  
  # Créer les plots de richness (sans log)
  richness_plots <- list()
  if(length(available_richness) > 0) {
    richness_plots <- create_richness_boxplots(data, trait_col, available_richness, group_colors)
  }
  
  # Combiner tous les plots
  all_plots <- c(emission_plots, richness_plots)
  
  # Sauvegarder les plots si demandé
  if(save_plots) {
    if(is.null(save_dir)) {
      save_dir <- here::here("figures", "graphs_screening")
    }
    
    # Créer le dossier s'il n'existe pas
    if(!dir.exists(save_dir)) {
      dir.create(save_dir, recursive = TRUE)
    }
    
    # Sauvegarder chaque plot
    for(plot_name in names(all_plots)) {
      # Nettoyer le nom du fichier
      file_name <- gsub("[^a-zA-Z0-9]", "_", plot_name)
      file_path <- file.path(save_dir, paste0(file_name, ".png"))
      
      # Sauvegarder avec une bonne résolution
      ggplot2::ggsave(
        filename = file_path,
        plot = all_plots[[plot_name]],
        width = 8,
        height = 6,
        dpi = 300,
        bg = "white"
      )
      
      message("Plot sauvegardé : ", file_path)
    }
  }
  
  return(data)
}

# Sous-fonction pour les émissions (avec log) - version simple
create_emission_boxplots<- function(data, trait_col, emission_cols, group_colors) {
  
  plots <- list()
  
  for(col in emission_cols) {
    
    plot_data <- data %>%
      dplyr::select(all_of(c(trait_col, col))) %>%
      dplyr::filter(!is.na(!!sym(col)))
    
    # Appliquer log10(x + 1)
    plot_data <- plot_data %>%
      dplyr::mutate(y_log = log10(!!sym(col) + 1))
    
    # Labels simples
    y_labels <- c(
      "Isoprene" = "log10(Isoprene + 1)",
      "Monoterpenes" = "log10(Monoterpenes + 1)",
      "Oxygenated-monoterpenes" = "log10(Oxygenated monoterpenes + 1)",
      "Total_monoterpenes" = "log10(Total monoterpenes + 1)"
    )
    
    # Créer le plot
    p <- ggplot2::ggplot(
      plot_data, 
      ggplot2::aes(x = !!sym(trait_col), y = y_log, fill = !!sym(trait_col))
    ) +
      ggplot2::geom_boxplot(
        alpha = 0.8,
        outlier.shape = NA,
        color = "black",
        linewidth = 0.5
      ) +
      ggplot2::geom_jitter(
        width = 0.15,
        height = 0,
        alpha = 0.6,
        size = 2,
        shape = 21,
        fill = "white",
        color = "black"
      ) +
      ggplot2::scale_fill_manual(values = group_colors) +
      ggplot2::labs(x = "", y = y_labels[col]) +
      ggplot2::theme_minimal() +
      ggplot2::theme(
        axis.text.y = ggplot2::element_text(size = 10),
        axis.title.y = ggplot2::element_text(size = 11, margin = ggplot2::margin(r = 10)),
        axis.text.x = ggplot2::element_text(size = 11, face = "bold"),
        plot.title = ggplot2::element_blank(),
        legend.position = "none",
        panel.grid.major.x = ggplot2::element_blank(),
        panel.grid.minor.x = ggplot2::element_blank(),
        panel.grid.major.y = ggplot2::element_line(color = "gray90", linewidth = 0.3),
        panel.border = ggplot2::element_rect(fill = NA, color = "gray70", linewidth = 0.5),
        plot.margin = ggplot2::margin(15, 15, 15, 15)
      ) +
      ggplot2::scale_y_continuous(
        breaks = scales::pretty_breaks(n = 5),
        labels = function(x) sprintf("%.2f", x),
        expand = ggplot2::expansion(mult = c(0.05, 0.1))
      )
    
    plots[[col]] <- p
  }
  
  return(plots)
}

# Sous-fonction pour les richness (sans log)
create_richness_boxplots <- function(data, trait_col, richness_cols, group_colors) {
  
  # Labels Y pour richness
  y_labels <- c(
    "Count_Monoterpenes" = "Monoterpenes richness",
    "Count_Oxygenated-monoterpenes" = "Oxygenated monoterpenes richness", 
    "Count_Total_monoterpenes" = "Total monoterpenes richness",
    "Count_Sesquiterpenes" = "Sesquiterpenes richness",
    "Count_Oxygenated-sesquiterpenes" = "Oxygenated sesquiterpenes richness"
  )
  
  plots <- list()
  
  for(col in richness_cols) {
    
    plot_data <- data %>%
      dplyr::select(all_of(c(trait_col, col))) %>%
      dplyr::filter(!is.na(!!sym(col)))
    
    # Créer le plot
    p <- ggplot2::ggplot(
      plot_data, 
      ggplot2::aes(x = !!sym(trait_col), y = !!sym(col), fill = !!sym(trait_col))
    ) +
      ggplot2::geom_boxplot(
        alpha = 0.8,
        outlier.shape = NA,
        color = "black",
        linewidth = 0.5
      ) +
      ggplot2::geom_jitter(
        width = 0.15,
        height = 0,
        alpha = 0.6,
        size = 2,
        shape = 21,
        fill = "white",
        color = "black"
      ) +
      ggplot2::scale_fill_manual(values = group_colors) +
      ggplot2::labs(x = "", y = y_labels[col]) +
      ggplot2::theme_minimal() +
      ggplot2::theme(
        axis.text.y = ggplot2::element_text(size = 10),
        axis.title.y = ggplot2::element_text(size = 11, margin = ggplot2::margin(r = 10)),
        axis.text.x = ggplot2::element_text(size = 11, face = "bold"),
        plot.title = ggplot2::element_blank(),
        legend.position = "none",
        panel.grid.major.x = ggplot2::element_blank(),
        panel.grid.minor.x = ggplot2::element_blank(),
        panel.grid.major.y = ggplot2::element_line(color = "gray90", linewidth = 0.3),
        panel.border = ggplot2::element_rect(fill = NA, color = "gray70", linewidth = 0.5),
        plot.margin = ggplot2::margin(15, 15, 15, 15)
      ) +
      ggplot2::scale_y_continuous(
        breaks = scales::pretty_breaks(),
        expand = ggplot2::expansion(mult = c(0.05, 0.1))
      )
    
    plots[[col]] <- p
  }
  
  return(plots)
}


create_stacked_barchart <- function(compound_mean_spagg, spagg_code) {
  
  # 1. Filtrer pour garder seulement les classes d'intérêt
  filtered_data <- compound_mean_spagg %>%
    dplyr::filter(class %in% c("Isoprene", "Monoterpenes", "Oxygenated-monoterpenes")) %>%
    dplyr::select(compound, class, starts_with("mean_"))
  
  # Extraire la colonne pour l'espèce spécifique
  spagg_col <- paste0("mean_", spagg_code)
  
  if(!spagg_col %in% names(filtered_data)) {
    stop(paste("L'espèce", spagg_code, "n'est pas présente dans les données"))
  }
  
  # 2. Préparer les données pour cette espèce
  plot_data <- filtered_data %>%
    dplyr::select(compound, class, !!sym(spagg_col)) %>%
    dplyr::mutate(
      # Extraire la valeur moyenne et l'écart-type
      value_str = !!sym(spagg_col),
      # Séparer moyenne et sd (format: "0.083 (0.032)")
      mean_value = as.numeric(stringr::str_extract(value_str, "^[0-9.]+")),
      sd_value = as.numeric(stringr::str_extract(value_str, "(?<=\\()[0-9.]+(?=\\))")),
      # Remplacer NA par 0
      mean_value = ifelse(is.na(mean_value), 0, mean_value),
      sd_value = ifelse(is.na(sd_value), 0, sd_value)
    ) %>%
    dplyr::filter(mean_value > 0) %>%  # Garder seulement les compounds avec des valeurs > 0
    dplyr::arrange(class, desc(mean_value)) %>%  # Trier
    dplyr::mutate(
      compound = factor(compound, levels = compound),  # Garder l'ordre
      # Calculer la position cumulée pour les barres d'erreur
      ymin_cum = cumsum(mean_value) - mean_value,  # Position de départ
      ymax_cum = cumsum(mean_value) + sd_value     # Position de fin (moyenne + SD)
    )
  
  # Vérifier s'il y a des données
  if(nrow(plot_data) == 0) {
    warning(paste("Aucune donnée pour l'espèce", spagg_code))
    return(NULL)
  }
  
  # 3. Créer une palette de couleurs
  n_compounds <- nrow(plot_data)
  
  if(n_compounds <= 8) {
    colors <- RColorBrewer::brewer.pal(8, "Set2")
  } else if(n_compounds <= 12) {
    colors <- RColorBrewer::brewer.pal(12, "Set3")
  } else {
    colors <- colorRampPalette(RColorBrewer::brewer.pal(12, "Set3"))(n_compounds)
  }
  
  # Assigner les couleurs
  plot_data <- plot_data %>%
    dplyr::mutate(color = colors[1:n_compounds])
  
  # 4. Créer le plot empilé
  p <- ggplot2::ggplot(
    plot_data,
    ggplot2::aes(x = spagg_code, y = mean_value, fill = compound)
  ) +
    # Barres empilées
    ggplot2::geom_bar(
      stat = "identity",
      position = ggplot2::position_stack(reverse = TRUE),
      color = "black",
      linewidth = 0.3,
      width = 0.6
    ) +
    # Points pour les extrémités des T-bars (hauts)
    ggplot2::geom_point(
      ggplot2::aes(y = ymax_cum),
      position = ggplot2::position_stack(reverse = TRUE),
      shape = 95,  # Tiret bas
      size = 8,
      color = "black",
      stroke = 1.5
    ) +
    # Segments verticaux pour les T-bars
    ggplot2::geom_segment(
      ggplot2::aes(
        x = as.numeric(factor(spagg_code)) - 0.15,
        xend = as.numeric(factor(spagg_code)) + 0.15,
        y = ymax_cum,
        yend = ymax_cum
      ),
      position = ggplot2::position_stack(reverse = TRUE),
      color = "black",
      linewidth = 0.8
    ) +
    # Segments horizontaux pour les T-bars
    ggplot2::geom_segment(
      ggplot2::aes(
        x = as.numeric(factor(spagg_code)),
        xend = as.numeric(factor(spagg_code)),
        y = cumsum(mean_value),  # Sommet de la barre
        yend = ymax_cum          # Sommet + SD
      ),
      position = ggplot2::position_stack(reverse = TRUE),
      color = "black",
      linewidth = 0.5
    ) +
    # Échelle de couleurs
    ggplot2::scale_fill_manual(
      values = setNames(plot_data$color, plot_data$compound),
      name = "Compound"
    ) +
    # Labels
    ggplot2::labs(
      title = paste("Composition BVOC -", spagg_code),
      x = "",
      y = expression("Emission rate" ~ (µg~g^{-1}~h^{-1})),
      caption = "Barres T : ±1 SD"
    ) +
    ggplot2::theme_minimal() +
    ggplot2::theme(
      plot.title = ggplot2::element_text(
        hjust = 0.5,
        size = 14,
        face = "bold",
        margin = ggplot2::margin(b = 15)
      ),
      axis.title.y = ggplot2::element_text(
        size = 12,
        margin = ggplot2::margin(r = 10)
      ),
      axis.text.y = ggplot2::element_text(size = 10),
      axis.text.x = ggplot2::element_text(
        size = 12,
        face = "bold"
      ),
      legend.title = ggplot2::element_text(size = 11, face = "bold"),
      legend.text = ggplot2::element_text(size = 9),
      legend.key.size = ggplot2::unit(0.6, "cm"),
      panel.grid.major.x = ggplot2::element_blank(),
      panel.grid.minor.x = ggplot2::element_blank(),
      panel.grid.major.y = ggplot2::element_line(color = "gray90", linewidth = 0.3),
      panel.grid.minor.y = ggplot2::element_blank(),
      panel.border = ggplot2::element_rect(fill = NA, color = "gray70", linewidth = 0.5),
      plot.margin = ggplot2::margin(15, 15, 15, 15),
      plot.caption = ggplot2::element_text(size = 9, color = "gray50")
    ) +
    ggplot2::scale_y_continuous(
      expand = ggplot2::expansion(mult = c(0, 0.1))
    )
  
  return(p)
}


create_all_stacked_barcharts <- function(compound_mean_spagg, save_dir = NULL) {
  
  if(is.null(save_dir)) {
    save_dir <- here::here("figures", "graphs_screening", "compo_species")
  }
  
  # Créer le dossier
  if(!dir.exists(save_dir)) {
    dir.create(save_dir, recursive = TRUE, showWarnings = FALSE)
    cat("Dossier créé:", save_dir, "\n")
  }
  
  # Créer une palette de couleurs cohérente pour tous les compounds
  all_compounds <- compound_mean_spagg %>%
    dplyr::filter(class %in% c("Isoprene", "Monoterpenes", "Oxygenated-monoterpenes")) %>%
    dplyr::pull(compound) %>%
    unique() %>%
    sort()
  
  n_compounds <- length(all_compounds)
  
  if(n_compounds <= 12) {
    base_colors <- RColorBrewer::brewer.pal(12, "Set3")
  } else {
    base_colors <- colorRampPalette(RColorBrewer::brewer.pal(12, "Set3"))(n_compounds)
  }
  
  color_palette <- setNames(base_colors[1:n_compounds], all_compounds)
  cat("Palette créée pour", n_compounds, "compounds\n")
  
  # Extraire les espèces
  spagg_cols <- names(compound_mean_spagg)[grepl("^mean_", names(compound_mean_spagg))]
  spagg_codes <- gsub("^mean_", "", spagg_cols)
  
  plots <- list()
  
  for(spagg in spagg_codes) {
    cat("Traitement de", spagg, "... ")
    
    # Utiliser la même palette de couleurs pour toutes les espèces
    p <- create_stacked_barchart(compound_mean_spagg, spagg, color_palette)
    
    if(!is.null(p)) {
      # Sauvegarder
      file_path <- file.path(save_dir, paste0("composition_", spagg, ".png"))
      
      ggplot2::ggsave(
        filename = file_path,
        plot = p,
        width = 12,  # Plus large pour la légende
        height = 8,
        dpi = 300,
        bg = "white"
      )
      
      plots[[spagg]] <- p
      cat("✓ Sauvegardé\n")
    } else {
      cat("✗ Pas de données\n")
    }
  }
  
  # Sauvegarder aussi la légende séparément
  if(length(plots) > 0) {
    save_legend_separately(color_palette, save_dir)
  }
  
  return(list(plots = plots, color_palette = color_palette))
}

# Fonction pour sauvegarder la légende séparément
save_legend_separately <- function(color_palette, save_dir) {
  # Créer un petit plot juste pour la légende
  legend_data <- data.frame(
    compound = names(color_palette),
    value = 1,
    stringsAsFactors = FALSE
  )
  
  p_legend <- ggplot2::ggplot(
    legend_data,
    ggplot2::aes(x = compound, y = value, fill = compound)
  ) +
    ggplot2::geom_bar(stat = "identity") +
    ggplot2::scale_fill_manual(values = color_palette, name = "Compounds") +
    ggplot2::theme(
      legend.position = "right",
      legend.box = "vertical",
      legend.key.size = ggplot2::unit(0.8, "cm"),
      legend.text = ggplot2::element_text(size = 9)
    ) +
    ggplot2::guides(fill = ggplot2::guide_legend(ncol = 2))
  
  # Extraire la légende
  legend <- cowplot::get_legend(p_legend)
  
  # Sauvegarder
  legend_file <- file.path(save_dir, "color_legend.png")
  ggplot2::ggsave(
    filename = legend_file,
    plot = legend,
    width = 6,
    height = max(8, length(color_palette) * 0.3),
    dpi = 300,
    bg = "white"
  )
  
  cat("Légende sauvegardée:", legend_file, "\n")
}



create_stacked_barchart_errorbar <- function(compound_mean_spagg, spagg_code, color_palette = NULL) {
  
  # Même préparation des données
  filtered_data <- compound_mean_spagg %>%
    dplyr::filter(class %in% c("Isoprene", "Monoterpenes", "Oxygenated-monoterpenes")) %>%
    dplyr::select(compound, class, starts_with("mean_"))
  
  spagg_col <- paste0("mean_", spagg_code)
  
  if(!spagg_col %in% names(filtered_data)) {
    warning(paste("L'espèce", spagg_code, "n'est pas présente"))
    return(NULL)
  }
  
  plot_data <- filtered_data %>%
    dplyr::select(compound, class, !!sym(spagg_col)) %>%
    dplyr::mutate(
      value_str = !!sym(spagg_col),
      mean_value = as.numeric(stringr::str_extract(value_str, "^[0-9.]+")),
      sd_value = as.numeric(stringr::str_extract(value_str, "(?<=\\()[0-9.]+(?=\\))"))
    ) %>%
    dplyr::mutate(
      mean_value = ifelse(is.na(mean_value), 0, mean_value),
      sd_value = ifelse(is.na(sd_value), 0, sd_value)
    ) %>%
    dplyr::filter(mean_value > 0) %>%
    dplyr::arrange(class, desc(mean_value))
  
  if(nrow(plot_data) == 0) {
    warning(paste("Aucune donnée pour", spagg_code))
    return(NULL)
  }
  
  # Calculer les positions cumulées pour geom_errorbar
  plot_data <- plot_data %>%
    dplyr::mutate(
      y_cum = cumsum(mean_value),  # Position cumulée (haut du segment)
      ymin_error = y_cum - sd_value,  # Bas de la barre d'erreur
      ymax_error = y_cum + sd_value   # Haut de la barre d'erreur
    )
  
  # Palette de couleurs cohérente
  all_compounds <- compound_mean_spagg %>%
    dplyr::filter(class %in% c("Isoprene", "Monoterpenes", "Oxygenated-monoterpenes")) %>%
    dplyr::pull(compound) %>%
    unique() %>%
    sort()
  
  n_compounds <- length(all_compounds)
  
  if(is.null(color_palette)) {
    if(n_compounds <= 12) {
      base_colors <- RColorBrewer::brewer.pal(12, "Set3")
    } else {
      base_colors <- colorRampPalette(RColorBrewer::brewer.pal(12, "Set3"))(n_compounds)
    }
    color_palette <- setNames(base_colors[1:n_compounds], all_compounds)
  }
  
  plot_colors <- color_palette[names(color_palette) %in% plot_data$compound]
  
  # Créer le plot avec geom_errorbar
  p <- ggplot2::ggplot(
    plot_data,
    ggplot2::aes(x = spagg_code, y = mean_value, fill = compound)
  ) +
    # Barres empilées
    ggplot2::geom_bar(
      stat = "identity",
      position = ggplot2::position_stack(reverse = FALSE),
      color = "black",
      linewidth = 0.3,
      width = 0.7
    ) +
    # Barres d'erreur T - positionnées au sommet de chaque segment
    ggplot2::geom_errorbar(
      ggplot2::aes(
        ymin = ymin_error,
        ymax = ymax_error
      ),
      position = ggplot2::position_stack(reverse = FALSE, vjust = 0),
      width = 0.2,
      linewidth = 0.8,
      color = "black"
    ) +
    ggplot2::scale_fill_manual(
      values = plot_colors,
      name = "Compound",
      breaks = all_compounds
    ) +
    ggplot2::labs(
      title = paste("Composition BVOC -", spagg_code),
      x = "",
      y = expression("Emission rate" ~ (µg~g^{-1}~h^{-1}))
    ) +
    ggplot2::theme_minimal() +
    ggplot2::theme(
      plot.title = ggplot2::element_text(hjust = 0.5, size = 14, face = "bold"),
      axis.text.x = ggplot2::element_text(size = 12, face = "bold")
    )
  
  return(p)
}



