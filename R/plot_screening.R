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
  
  # Vérifier que les colonnes nécessaires existent pour les sesquiterpènes
  # Si elles n'existent pas, les créer avec 0
  if("Sesquiterpenes" %in% names(df_combined) && "Oxygenated-sesquiterpenes" %in% names(df_combined)) {
    df_combined <- df_combined |> 
      dplyr::mutate(
        Total_sesquiterpenes = Sesquiterpenes + `Oxygenated-sesquiterpenes`
      )
  } else if("Sesquiterpenes" %in% names(df_combined)) {
    df_combined <- df_combined |> 
      dplyr::mutate(
        `Oxygenated-sesquiterpenes` = 0,
        Total_sesquiterpenes = Sesquiterpenes
      )
  } else if("Oxygenated-sesquiterpenes" %in% names(df_combined)) {
    df_combined <- df_combined |> 
      dplyr::mutate(
        Sesquiterpenes = 0,
        Total_sesquiterpenes = `Oxygenated-sesquiterpenes`
      )
  } else {
    df_combined <- df_combined |> 
      dplyr::mutate(
        Sesquiterpenes = 0,
        `Oxygenated-sesquiterpenes` = 0,
        Total_sesquiterpenes = 0
      )
  }
  
  # Calculer le compte total des monoterpènes et sesquiterpènes
  df_combined <- df_combined |> 
    dplyr::mutate(
      Count_Total_monoterpenes = Count_Monoterpenes + `Count_Oxygenated-monoterpenes`,
      Count_Total_sesquiterpenes = Count_Sesquiterpenes + `Count_Oxygenated-sesquiterpenes`
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
      Sesquiterpenes,
      `Oxygenated-sesquiterpenes`,
      Total_sesquiterpenes,
      Count_Monoterpenes,
      `Count_Oxygenated-monoterpenes`,
      Count_Total_monoterpenes,
      Count_Sesquiterpenes,
      `Count_Oxygenated-sesquiterpenes`,
      Count_Total_sesquiterpenes,
      Taxon, 
      PAR_algo, 
      T_algo_K,
      Standardized
    ) |> 
    dplyr::mutate(Taxon = stringr::str_replace_all(Taxon, " ", "_"))
  
  return(bvocs_samples_wide)
}



create_all_boxplots <- function(data, trait_col, 
                                remove_zeros = TRUE,
                                test_emissions = TRUE,
                                test_richness = FALSE,  # Nouveau paramètre
                                apply_log_transform = TRUE,
                                save_path = NULL) {
  df <- data
  
  # Setup
  data <- data |> 
    dplyr::filter(tolower(!!sym(trait_col)) != "intermediate")
  
  unique_groups <- unique(data[[trait_col]])
  if(length(unique_groups) > 2) {
    unique_groups <- unique_groups[1:2]
    data <- data |> 
      dplyr::filter(!!sym(trait_col) %in% unique_groups)
  }
  
  colors <- setNames(c("#bd1b6e", "#ff862f")[1:length(unique_groups)], unique_groups)
  
  # Create emission plots
  emission_results <- create_emission_boxplots(
    data = data,
    trait_col = trait_col,
    emission_cols = intersect(
      c("Isoprene", "Monoterpenes", "Oxygenated-monoterpenes", 
        "Total_monoterpenes", "Sesquiterpenes", "Oxygenated-sesquiterpenes",
        "Total_sesquiterpenes"),
      names(data)
    ),
    group_colors = colors,
    remove_zeros = remove_zeros,
    include_stat_test = test_emissions,
    apply_log_transform = apply_log_transform
  )
  
  emission_plots <- emission_results$plots
  
  # Create richness plots avec test optionnel
  richness_results <- create_richness_boxplots(
    data = data,
    trait_col = trait_col,
    richness_cols = intersect(
      c("Count_Monoterpenes", "Count_Oxygenated-monoterpenes", 
        "Count_Total_monoterpenes", "Count_Total_sesquiterpenes",
        "Count_Sesquiterpenes", "Count_Oxygenated-sesquiterpenes"),
      names(data)
    ),
    group_colors = colors,
    include_stat_test = test_richness
  )
  
  richness_plots <- richness_results$plots
  
  # Combine all plots
  all_plots <- c(emission_plots, richness_plots)
  
  # Combine all statistics
  all_statistics <- list()
  if(!is.null(emission_results$statistics)) {
    all_statistics <- c(all_statistics, emission_results$statistics)
  }
  if(!is.null(richness_results$statistics)) {
    all_statistics <- c(all_statistics, richness_results$statistics)
  }
  
  # Save plots if requested
  if(!is.null(save_path)) {
    if(!dir.exists(save_path)) {
      dir.create(save_path, recursive = TRUE, showWarnings = FALSE)
    }
    
    for(name in names(all_plots)) {
      file_name <- gsub("[^a-zA-Z0-9]", "_", name)
      file_path <- file.path(save_path, paste0(file_name, ".png"))
      
      ggplot2::ggsave(
        filename = file_path,
        plot = all_plots[[name]],
        width = 8,
        height = 6,
        dpi = 300,
        bg = "white"
      )
    }
    
    # Save all statistical results
    if(length(all_statistics) > 0) {
      stats_df <- do.call(rbind, lapply(all_statistics, function(x) {
        data.frame(
          variable = x$compound,
          groups = x$groups,
          p_value = x$p_value,
          p_value_display = ifelse(!is.null(x$p_value_display), x$p_value_display, NA),
          n_group1 = x$n_group1,
          n_group2 = x$n_group2,
          is_richness = grepl("Count_", x$compound),
          log_transformed = ifelse(grepl("Count_", x$compound), FALSE, apply_log_transform),
          stringsAsFactors = FALSE
        )
      }))
      
      write.csv(stats_df, file.path(save_path, "statistical_results.csv"), row.names = FALSE)
    }
  }
  
  return(df)

}



# Function to create 4 types of distributions for each compound in a line format
plot_emission_distributions <- function(data, emission_cols, 
                                          save_plots = FALSE, output_dir = NULL) {
  
  # Prepare for saving histograms if requested
  if(save_plots && !is.null(output_dir)) {
    if(!dir.exists(output_dir)) {
      dir.create(output_dir, recursive = TRUE)
    }
  }
  
  # List to store all individual plots
  all_plots <- list()
  
  # List to store line plots (4 plots per compound in one line)
  line_plots <- list()
  
  # For each compound, create the 4 distributions
  for(col in emission_cols) {
    
    # Prepare base data
    plot_data <- data %>%
      dplyr::select(all_of(col)) %>%
      dplyr::filter(!is.na(!!sym(col)))
    
    # Data without zeros
    plot_data_no_zeros <- plot_data %>%
      dplyr::filter(!!sym(col) > 0)
    
    # Calculate zero statistics
    total_count <- nrow(plot_data)
    zero_count <- sum(plot_data[[col]] == 0, na.rm = TRUE)
    zero_percentage <- ifelse(total_count > 0, round(zero_count / total_count * 100, 1), 0)
    
    # Internal function to create a histogram with consistent styling
    create_histogram <- function(data_vec, title, subtitle = "", x_label = "Value") {
      
      # Create dataframe for ggplot
      hist_data <- data.frame(value = data_vec)
      
      # Create the plot
      p <- ggplot2::ggplot(hist_data, ggplot2::aes(x = value)) +
        ggplot2::geom_histogram(
          fill = "#6baed6",
          color = "black",
          alpha = 0.7,
          bins = 30
        ) +
        ggplot2::labs(
          x = x_label,
          y = "Frequency",
          title = title,
          subtitle = subtitle
        ) +
        ggplot2::theme_minimal() +
        ggplot2::theme(
          plot.title = ggplot2::element_text(size = 9, face = "bold", hjust = 0.5),
          plot.subtitle = ggplot2::element_text(size = 7, hjust = 0.5, color = "gray50"),
          panel.grid.minor = ggplot2::element_blank(),
          panel.border = ggplot2::element_rect(fill = NA, color = "gray70", linewidth = 0.5),
          axis.text = ggplot2::element_text(size = 6),
          axis.title = ggplot2::element_text(size = 7)
        )
      
      return(p)
    }
    
    # Create the 4 plots
    plots_list <- list()
    
    # 1. Raw distribution (with zeros)
    if(nrow(plot_data) > 0) {
      plots_list[["raw_with_zeros"]] <- create_histogram(
        data_vec = plot_data[[col]],
        title = "Raw (with zeros)",
        subtitle = paste("n =", total_count),
        x_label = "Value"
      )
    } else {
      plots_list[["raw_with_zeros"]] <- ggplot2::ggplot() + 
        ggplot2::annotate("text", x = 0.5, y = 0.5, label = "No data", size = 3) +
        ggplot2::theme_void()
    }
    
    # 2. Log transformed distribution (with zeros)
    if(nrow(plot_data) > 0) {
      log_values <- log10(plot_data[[col]] + 1)
      plots_list[["log_with_zeros"]] <- create_histogram(
        data_vec = log_values,
        title = "Log(x+1) (with zeros)",
        subtitle = paste("n =", total_count),
        x_label = "log10(x+1)"
      )
    } else {
      plots_list[["log_with_zeros"]] <- ggplot2::ggplot() + 
        ggplot2::annotate("text", x = 0.5, y = 0.5, label = "No data", size = 3) +
        ggplot2::theme_void()
    }
    
    # 3. Raw distribution without zeros
    if(nrow(plot_data_no_zeros) > 0) {
      plots_list[["raw_no_zeros"]] <- create_histogram(
        data_vec = plot_data_no_zeros[[col]],
        title = "Raw (no zeros)",
        subtitle = paste("n =", nrow(plot_data_no_zeros)),
        x_label = "Value"
      )
    } else {
      plots_list[["raw_no_zeros"]] <- ggplot2::ggplot() + 
        ggplot2::annotate("text", x = 0.5, y = 0.5, label = "No non-zero data", size = 3) +
        ggplot2::theme_void()
    }
    
    # 4. Log transformed distribution without zeros
    if(nrow(plot_data_no_zeros) > 0) {
      log_values_no_zeros <- log10(plot_data_no_zeros[[col]] + 1)
      plots_list[["log_no_zeros"]] <- create_histogram(
        data_vec = log_values_no_zeros,
        title = "Log(x+1) (no zeros)",
        subtitle = paste("n =", nrow(plot_data_no_zeros)),
        x_label = "log10(x+1)"
      )
    } else {
      plots_list[["log_no_zeros"]] <- ggplot2::ggplot() + 
        ggplot2::annotate("text", x = 0.5, y = 0.5, label = "No non-zero data", size = 3) +
        ggplot2::theme_void()
    }
    
    # Combine the 4 plots into one line
    line_plot <- patchwork::wrap_plots(
      plots_list[["raw_with_zeros"]],
      plots_list[["log_with_zeros"]],
      plots_list[["raw_no_zeros"]],
      plots_list[["log_no_zeros"]],
      ncol = 4,
      nrow = 1
    ) +
      patchwork::plot_annotation(
        title = paste("Distribution Analysis -", col),
        subtitle = paste("Total samples:", total_count, "| Zeros:", zero_count, 
                         "(", zero_percentage, "%)"),
        theme = ggplot2::theme(
          plot.title = ggplot2::element_text(size = 12, face = "bold", hjust = 0.5),
          plot.subtitle = ggplot2::element_text(size = 9, hjust = 0.5)
        )
      )
    
    # Store individual plots
    all_plots[[paste0(col, "_raw_with_zeros")]] <- plots_list[["raw_with_zeros"]]
    all_plots[[paste0(col, "_log_with_zeros")]] <- plots_list[["log_with_zeros"]]
    all_plots[[paste0(col, "_raw_no_zeros")]] <- plots_list[["raw_no_zeros"]]
    all_plots[[paste0(col, "_log_no_zeros")]] <- plots_list[["log_no_zeros"]]
    
    # Store line plot
    line_plots[[col]] <- line_plot
    
    # Save line plot if requested
    if(save_plots && !is.null(output_dir)) {
      ggplot2::ggsave(
        filename = file.path(output_dir, paste0(col, "_distribution_line.png")),
        plot = line_plot,
        width = 12,  # Wider for horizontal layout
        height = 4,  # Shorter height
        dpi = 300
      )
    }
  }
  
  # Create A4-sized combined plot with all compounds
  if(length(emission_cols) > 0) {
    # Calculate A4 dimensions (in inches)
    a4_width_in <- 8.27  # A4 width in inches
    a4_height_in <- 11.69  # A4 height in inches
    
    # Calculate space needed for each compound
    # Allow more height per compound if we have fewer compounds
    n_compounds <- length(emission_cols)
    compound_height <- min(2.5, a4_height_in / n_compounds)
    
    # Combine all line plots vertically
    combined_a4_plot <- patchwork::wrap_plots(
      line_plots,
      ncol = 1,
      heights = rep(compound_height, n_compounds)
    ) +
      patchwork::plot_annotation(
        title = "Compound Emission Distributions",
        subtitle = paste("Analysis of", length(emission_cols), "compounds"),
        theme = ggplot2::theme(
          plot.title = ggplot2::element_text(size = 14, face = "bold", hjust = 0.5, 
                                             margin = ggplot2::margin(b = 10)),
          plot.subtitle = ggplot2::element_text(size = 11, hjust = 0.5, 
                                                margin = ggplot2::margin(b = 15))
        )
      )
    
    # Save A4 combined plot if requested
    if(save_plots && !is.null(output_dir)) {
      ggplot2::ggsave(
        filename = file.path(output_dir, "all_compounds_A4_distributions.pdf"),
        plot = combined_a4_plot,
        width = a4_width_in,
        height = a4_height_in,
        device = "pdf"
      )
      
      # Also save as PNG
      ggplot2::ggsave(
        filename = file.path(output_dir, "all_compounds_A4_distributions.png"),
        plot = combined_a4_plot,
        width = a4_width_in,
        height = a4_height_in,
        dpi = 300
      )
    }
    
    # Store A4 plot
    all_plots[["A4_combined_plot"]] <- combined_a4_plot
  }
  
  if(save_plots && !is.null(output_dir)) {
    cat("Distribution plots saved in:", output_dir, "\n")
  }
  
  # Return both individual and combined plots
  return(list(
    individual_plots = all_plots,
    line_plots = line_plots,
    A4_combined_plot = if(exists("combined_a4_plot")) combined_a4_plot else NULL
  ))
}




# Function to create emission boxplots with statistical tests
create_emission_boxplots <- function(data, trait_col, emission_cols, group_colors,
                                     remove_zeros = FALSE, 
                                     include_stat_test = FALSE,
                                     apply_log_transform = TRUE) {
  
  plots <- list()
  stat_results <- list()
  
  for(col in emission_cols) {
    
    # 1. Préparer les données
    plot_data <- data %>%
      dplyr::select(all_of(c(trait_col, col))) %>%
      dplyr::filter(!is.na(!!sym(col)))
    
    if(remove_zeros) {
      plot_data <- plot_data %>% dplyr::filter(!!sym(col) > 0)
    }
    
    # 2. Appliquer transformation log SI demandée
    if(apply_log_transform) {
      plot_data <- plot_data %>%
        dplyr::mutate(y_value = log10(!!sym(col) + 1))
    } else {
      plot_data <- plot_data %>%
        dplyr::mutate(y_value = !!sym(col))
    }
    
    # 3. Obtenir les groupes
    groups <- unique(plot_data[[trait_col]])
    
    # 4. Calculer statistiques des zéros
    total_original <- nrow(data %>% dplyr::filter(!is.na(!!sym(col))))
    zero_count <- sum(data[[col]] == 0, na.rm = TRUE)
    zero_percentage <- ifelse(total_original > 0, round(zero_count / total_original * 100, 1), 0)
    
    # 5. Créer le titre
    title_parts <- c(col)
    
    if(apply_log_transform) {
      title_parts <- c(title_parts, "(log10 scale)")
    }
    
    title_text <- paste(title_parts, collapse = " ")
    
    # 6. Sous-titre avec infos échantillon
    subtitle_parts <- c(paste("n =", nrow(plot_data)))
    
    if(!remove_zeros && zero_count > 0) {
      subtitle_parts <- c(subtitle_parts, paste(zero_count, "zeros"))
    }
    
    if(remove_zeros) {
      subtitle_parts <- c(subtitle_parts, "zeros excluded")
    }
    
    subtitle_text <- paste(subtitle_parts, collapse = " | ")
    
    # 7. Test statistique
    p_value_display <- NULL
    
    if(include_stat_test && length(groups) == 2) {
      
      # Utiliser les mêmes données que pour le plot
      group1_data <- plot_data %>% 
        dplyr::filter(!!sym(trait_col) == groups[1]) %>% 
        dplyr::pull(y_value)
      
      group2_data <- plot_data %>% 
        dplyr::filter(!!sym(trait_col) == groups[2]) %>% 
        dplyr::pull(y_value)
      
      if(length(group1_data) >= 3 && length(group2_data) >= 3) {
        test_result <- wilcox.test(group1_data, group2_data, exact = FALSE)
        p_value <- test_result$p.value
        
        # Formater la p-valeur (simplifiée)
        if(p_value < 0.001) {
          p_value_display <- "p < 0.001"
        } else {
          p_value_display <- paste("p =", format(round(p_value, 3), nsmall = 3))
        }
        
        stat_results[[col]] <- list(
          compound = col,
          groups = paste(groups[1], "vs", groups[2]),
          p_value = p_value,
          p_value_display = p_value_display,
          n_group1 = length(group1_data),
          n_group2 = length(group2_data)
        )
      }
    }
    
    # 8. Définir le label y selon le composé
    if(col == "Isoprene") {
      if(apply_log_transform) {
        y_label <- bquote(log[10] ~ ("Isoprene EF" ~ (µg~g^{-1}~h^{-1}) + 1))
      } else {
        y_label <- expression("Isoprene EF" ~ (µg~g^{-1}~h^{-1}))
      }
    } else if(col == "Monoterpenes") {
      if(apply_log_transform) {
        y_label <- bquote(log[10] ~ ("Monoterpenes EF" ~ (µg~g^{-1}~h^{-1}) + 1))
      } else {
        y_label <- expression("Monoterpenes EF" ~ (µg~g^{-1}~h^{-1}))
      }
    } else if(col == "Total_monoterpenes") {
      if(apply_log_transform) {
        y_label <- bquote(log[10] ~ ("Monoterpenes EF" ~ (µg~g^{-1}~h^{-1}) + 1))
      } else {
        y_label <- expression("Monoterpenes EF" ~ (µg~g^{-1}~h^{-1}))
      }
    } else if(col == "Sesquiterpenes") {
      if(apply_log_transform) {
        y_label <- bquote(log[10] ~ ("Sesquiterpenes EF" ~ (µg~g^{-1}~h^{-1}) + 1))
      } else {
        y_label <- expression("Sesquiterpenes EF" ~ (µg~g^{-1}~h^{-1}))
      }
    } else if(col == "Total_sesquiterpenes") {
      if(apply_log_transform) {
        y_label <- bquote(log[10] ~ ("Sesquiterpenes EF" ~ (µg~g^{-1}~h^{-1}) + 1))
      } else {
        y_label <- expression("Sesquiterpenes EF" ~ (µg~g^{-1}~h^{-1}))
      }
    } else if(col == "Oxygenated-monoterpenes") {
      if(apply_log_transform) {
        y_label <- bquote(log[10] ~ ("Oxygenated monoterpenes EF" ~ (µg~g^{-1}~h^{-1}) + 1))
      } else {
        y_label <- expression("Oxygenated monoterpenes EF" ~ (µg~g^{-1}~h^{-1}))
      }
    } else if(col == "Oxygenated-sesquiterpenes") {
      if(apply_log_transform) {
        y_label <- bquote(log[10] ~ ("Oxygenated sesquiterpenes EF" ~ (µg~g^{-1}~h^{-1}) + 1))
      } else {
        y_label <- expression("Oxygenated sesquiterpenes EF" ~ (µg~g^{-1}~h^{-1}))
      }
    } else {
      # Fallback pour les autres colonnes
      if(apply_log_transform) {
        y_label <- paste0("log10(", col, " + 1)")
      } else {
        y_label <- col
      }
    }
    
    # 9. Créer le plot avec aes() correctement préfixé
    p <- ggplot2::ggplot(
      plot_data,
      ggplot2::aes(x = !!sym(trait_col), y = y_value, fill = !!sym(trait_col))
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
      ggplot2::labs(
        x = "",
        y = y_label,
        title = title_text,
        subtitle = subtitle_text
      ) +
      ggplot2::theme_minimal() +
      ggplot2::theme(
        axis.text.y = ggplot2::element_text(size = 10),
        axis.title.y = ggplot2::element_text(size = 11, margin = ggplot2::margin(r = 10)),
        axis.text.x = ggplot2::element_text(size = 11, face = "bold"),
        plot.title = ggplot2::element_text(size = 12, face = "bold", hjust = 0.5),
        plot.subtitle = ggplot2::element_text(size = 9, hjust = 0.5, color = "gray50"),
        legend.position = "none",
        panel.grid.major.x = ggplot2::element_blank(),
        panel.grid.minor.x = ggplot2::element_blank(),
        panel.grid.major.y = ggplot2::element_line(color = "gray90", linewidth = 0.3),
        panel.border = ggplot2::element_rect(fill = NA, color = "gray70", linewidth = 0.5),
        plot.margin = ggplot2::margin(15, 15, 15, 15)
      ) +
      ggplot2::scale_y_continuous(
        breaks = scales::pretty_breaks(n = 5)
      )
    
    # 10. Ajouter la p-valeur en noir si disponible
    if(!is.null(p_value_display)) {
      # Trouver la position y maximale
      y_max <- max(plot_data$y_value, na.rm = TRUE)
      y_min <- min(plot_data$y_value, na.rm = TRUE)
      y_range <- y_max - y_min
      
      # Position pour l'annotation (au-dessus du boxplot)
      annotation_y <- y_max + (0.08 * y_range)
      
      p <- p +
        ggplot2::annotate(
          "text",
          x = 1.5,  # Centre entre les deux groupes
          y = annotation_y,
          label = p_value_display,
          size = 4,
          color = "black",  # Noir
          fontface = "bold",
          vjust = 0
        )
    }
    
    plots[[col]] <- p
  }
  
  return(list(
    plots = plots,
    statistics = if(length(stat_results) > 0) stat_results else NULL
  ))
}


# Sous-fonction pour les richness (sans log)
create_richness_boxplots <- function(data, trait_col, richness_cols, group_colors,
                                        include_stat_test = FALSE,
                                        remove_zeros = FALSE) {
  
  # Labels pour richness avec expressions si besoin
  y_labels <- list(
    "Count_Monoterpenes" = "Monoterpenes richness",
    "Count_Oxygenated-monoterpenes" = "Oxygenated monoterpenes richness", 
    "Count_Total_monoterpenes" = "Monoterpenes richness",
    "Count_Sesquiterpenes" = "Sesquiterpenes richness",
    "Count_Oxygenated-sesquiterpenes" = "Oxygenated sesquiterpenes richness",
    "Count_Total_sesquiterpenes" = "Sesquiterpenes richness"
  )
  
  plots <- list()
  stat_results <- list()
  
  for(col in richness_cols) {
    
    # 1. Préparer les données
    plot_data <- data %>%
      dplyr::select(all_of(c(trait_col, col))) %>%
      dplyr::filter(!is.na(!!sym(col)))
    
    # Option: remove zeros
    if(remove_zeros) {
      plot_data <- plot_data %>% dplyr::filter(!!sym(col) > 0)
    }
    
    # 2. Utiliser les données directement
    plot_data <- plot_data %>%
      dplyr::mutate(y_value = !!sym(col))
    
    # 3. Obtenir les groupes
    groups <- unique(plot_data[[trait_col]])
    
    # 4. Calculer statistiques
    total_original <- nrow(data %>% dplyr::filter(!is.na(!!sym(col))))
    zero_count <- sum(data[[col]] == 0, na.rm = TRUE)
    zero_percentage <- ifelse(total_original > 0, round(zero_count / total_original * 100, 1), 0)
    
    # 5. Créer le titre
    title_text <- gsub("Count_", "", col)
    
    # 6. Sous-titre avec infos échantillon
    subtitle_parts <- c(paste("n =", nrow(plot_data)))
    
    if(!remove_zeros && zero_count > 0) {
      subtitle_parts <- c(subtitle_parts, paste(zero_count, "zeros"))
    }
    
    if(remove_zeros) {
      subtitle_parts <- c(subtitle_parts, "zeros excluded")
    }
    
    subtitle_text <- paste(subtitle_parts, collapse = " | ")
    
    # 7. Test statistique
    p_value_display <- NULL
    
    if(include_stat_test && length(groups) == 2) {
      
      # Utiliser les mêmes données que pour le plot
      group1_data <- plot_data %>% 
        dplyr::filter(!!sym(trait_col) == groups[1]) %>% 
        dplyr::pull(y_value)
      
      group2_data <- plot_data %>% 
        dplyr::filter(!!sym(trait_col) == groups[2]) %>% 
        dplyr::pull(y_value)
      
      if(length(group1_data) >= 3 && length(group2_data) >= 3) {
        test_result <- wilcox.test(group1_data, group2_data, exact = FALSE)
        p_value <- test_result$p.value
        
        # Formater la p-valeur
        if(p_value < 0.001) {
          p_value_display <- "p < 0.001"
        } else {
          p_value_display <- paste("p =", format(round(p_value, 3), nsmall = 3))
        }
        
        stat_results[[col]] <- list(
          compound = col,
          groups = paste(groups[1], "vs", groups[2]),
          p_value = p_value,
          p_value_display = p_value_display,
          n_group1 = length(group1_data),
          n_group2 = length(group2_data)
        )
      }
    }
    
    # 8. Définir le label y
    if(col %in% names(y_labels)) {
      y_label <- y_labels[[col]]
    } else {
      y_label <- gsub("Count_", "", col)
    }
    
    # 9. Créer le plot
    p <- ggplot2::ggplot(
      plot_data,
      ggplot2::aes(x = !!sym(trait_col), y = y_value, fill = !!sym(trait_col))
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
      ggplot2::labs(
        x = "",
        y = y_label,
        title = title_text,
        subtitle = subtitle_text
      ) +
      ggplot2::theme_minimal() +
      ggplot2::theme(
        axis.text.y = ggplot2::element_text(size = 10),
        axis.title.y = ggplot2::element_text(size = 11, margin = ggplot2::margin(r = 10)),
        axis.text.x = ggplot2::element_text(size = 11, face = "bold"),
        plot.title = ggplot2::element_text(size = 12, face = "bold", hjust = 0.5),
        plot.subtitle = ggplot2::element_text(size = 9, hjust = 0.5, color = "gray50"),
        legend.position = "none",
        panel.grid.major.x = ggplot2::element_blank(),
        panel.grid.minor.x = ggplot2::element_blank(),
        panel.grid.major.y = ggplot2::element_line(color = "gray90", linewidth = 0.3),
        panel.border = ggplot2::element_rect(fill = NA, color = "gray70", linewidth = 0.5),
        plot.margin = ggplot2::margin(15, 15, 15, 15)
      ) +
      ggplot2::scale_y_continuous(
        breaks = scales::pretty_breaks(n = 5),
        expand = ggplot2::expansion(mult = c(0.05, ifelse(!is.null(p_value_display), 0.12, 0.1)))
      )
    
    # 10. Ajouter la p-valeur en noir si disponible
    if(!is.null(p_value_display)) {
      y_max <- max(plot_data$y_value, na.rm = TRUE)
      
      p <- p +
        ggplot2::annotate(
          "text",
          x = 1.5,  # Centre entre les deux groupes
          y = y_max * 1.08,
          label = p_value_display,
          size = 4,
          color = "black",
          fontface = "bold",
          vjust = 0
        )
    }
    
    plots[[col]] <- p
  }
  
  return(list(
    plots = plots,
    statistics = if(length(stat_results) > 0) stat_results else NULL
  ))
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



