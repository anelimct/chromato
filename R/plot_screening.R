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




#' Title
#'
#' @param data 
#' @param trait_col 
#' @param remove_zeros 
#' @param test_emissions 
#' @param test_richness 
#' @param apply_log_transform 
#' @param save_path 
#'
#' @returns box-plots emission EF and richness according to leaf pheno decidious vs evergreen for the screening species
#' @export
#'
#' @examples
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





#' Title
#'
#' @param data 
#'
#' @returns data formated to do multivariate analysis for the screening visualization
#' @export
#'
#' @examples
prepare_plsda_data <- function(data) {
  # Step 1: Filter rows where class is NOT "Oxygenated-sesquiterpenes" or "Sesquiterpenes"
  filtered_data <- data |> 
    dplyr::filter(!class %in% c("Oxygenated-sesquiterpenes", "Sesquiterpenes"))
  
  # Step 2: Remove specified columns (class, smiles, inchikey)
  filtered_data <- filtered_data |> 
    dplyr::select(-class, -smiles, -inchikey)
  
  # Step 3: Set compound names as row names and remove the compound column
  # Then transpose the data
  # First convert to regular data frame if it's a tibble
  filtered_df <- as.data.frame(filtered_data)
  
  # Set row names to compound names and remove the column
  rownames(filtered_df) <- filtered_df$compound
  filtered_df$compound <- NULL
  
  # Step 4: Transpose the data (compounds become columns, samples become rows)
  transposed_data <- as.data.frame(t(filtered_df))
  
  # Convert all columns to numeric (they might be character due to NAs)
  transposed_data <- transposed_data |> 
    dplyr::mutate(across(everything(), as.numeric)) |> 
    dplyr::mutate(across(everything(), ~ ifelse(is.na(.), 0, .))) 
  
  transposed_data_N <- transposed_data |> normaliser_dataframe()
  
  gragg <- factor(stringr::str_to_upper(substr(rownames(transposed_data), 1, 4)))
  
  return(list(data = transposed_data_N, gragg = gragg, data_chemodiv = transposed_data))
}

plot_plsda <- function(data, group){
  
  
}


extract_value_sd <- function(x) {
  if (is.na(x) || x == "") return(c(NA, NA))
  parts <- strsplit(x, " \\(|\\)")[[1]]
  value <- as.numeric(parts[1])
  sd <- as.numeric(parts[2])
  c(value, sd)
}

# Function to plot emissions for a given class (isoprene or monoterpenes)
create_emission_plot <- function(df, class_name, y_axis_label = "Emission Rate") {
  # Filter for the class
  df_class <- df %>%
    filter(class %in% class_name)
  
  # Get the names of species columns (assuming columns 5 to end)
  species_cols <- names(df_class)[5:ncol(df_class)]
  
  # Reshape data for plotting
  df_long <- df_class %>%
    pivot_longer(cols = all_of(species_cols), names_to = "species", values_to = "value_sd") %>%
    mutate(value = sapply(value_sd, function(x) extract_value_sd(x)[1]),
           sd = sapply(value_sd, function(x) extract_value_sd(x)[2])) %>%
    group_by(species) %>%
    summarise(value = sum(value, na.rm = TRUE),
              sd = sqrt(sum(sd^2, na.rm = TRUE)))  # Sum of variances for combined SD
  
  # Create and return the plot
  p <- ggplot(df_long, aes(x = species, y = value)) +
    geom_bar(stat = "identity", fill = "grey", color = "black") +
    geom_errorbar(aes(ymin = value - sd, ymax = value + sd), width = 0.2) +
    labs(x = NULL, y = y_axis_label) +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))
  return(p)  # Explicitly return the plot
}

# Function to create a cow plot (stacked vertically)
plot_cow_emissions <- function(df) {
  # Create individual plots
  isoprene_plot <- create_emission_plot(df, "Isoprene", y_axis_label = expression("Isoprene EF" ~ (mu*g ~ "." ~g^{-1} ~ "." ~h^{-1})))+ theme(axis.title.x = element_blank())
  monoterpenes_plot <- create_emission_plot(df, c("Monoterpenes", "Oxygenated-monoterpenes"), y_axis_label = expression("Monoterpenes EF" ~ (mu*g ~ "." ~g^{-1} ~ "." ~h^{-1})))
  
  # Combine plots vertically with shared x-axis
  combined_plot <- isoprene_plot +
    monoterpenes_plot +
    plot_layout(ncol = 1)  # Stack vertically
  print(combined_plot)
}


pie_chart_screening <- function(pie_chart_emission_screening){
  
  
  pie_chart_emission_screening |> dplyr::filter(class != "Oxygenated-sesquiterpenes")
  pie_chart_emission_screening |> dplyr::filter(class != "Sesquiterpenes")
  species_cols <- grep("^mean_", names(pie_chart_emission_screening), value = TRUE)
  
  # Convert species columns to numeric and handle NA values
  pie_chart_numeric <- pie_chart_emission_screening  |> 
    dplyr::mutate(across(all_of(species_cols), ~ as.numeric(ifelse(. == "NA", NA, .))))
  
  # Now group by class and sum for each species, then pivot
  result <- pie_chart_numeric |> 
    dplyr::group_by(class) |> 
    # Sum across all species columns, removing NAs
    dplyr::summarise(across(all_of(species_cols), ~ sum(., na.rm = TRUE))) |> 
    # Pivot to long format
    tidyr::pivot_longer(
      cols = all_of(species_cols),
      names_to = "species",
      values_to = "total_emission"
    ) |> 
    # Remove "mean_" prefix from species names if desired
    dplyr::mutate(species = gsub("^mean_", "", species)) %>%
    # Pivot back to wide format with classes as columns (optional, based on your description)
    tidyr::pivot_wider(
      names_from = class,
      values_from = total_emission,
      values_fill = 0
    ) |> dplyr::mutate(monoterpenes = `Monoterpenes` + `Oxygenated-monoterpenes`) |> dplyr::mutate(isoprene = `Isoprene`) |> dplyr::mutate(Sum = `isoprene` + `monoterpenes`)|>     dplyr::mutate(type = dplyr::if_else(!is.na(Sum), 
                                             ifelse(isoprene >= 1 & monoterpenes > 0.2, "both", 
                                                    ifelse(monoterpenes > 0.2, "mono", 
                                                           ifelse(isoprene > 1, "iso", "NE"))), 
                                             NA_character_)) |> dplyr::mutate(sub_type = dplyr::case_when(
    # For iso type
    type == "iso" & isoprene < 10 ~ "low",
    type == "iso" & isoprene <= 30 ~ "medium",
    type == "iso" ~ "high",
    # For mono type
    type == "mono" & monoterpenes < 2 ~ "low",
    type == "mono" & monoterpenes <= 5.1 ~ "medium",
    type == "mono" ~ "high",
    # For both type - you need to decide thresholds for both
    type == "both" & isoprene < 10 & monoterpenes < 2 ~ "low-low",
    type == "both" & isoprene < 10 & monoterpenes <= 5.1 ~ "low-medium",
    type == "both" & isoprene < 10 ~ "low-high",
    type == "both" & isoprene <= 30 & monoterpenes < 2 ~ "medium-low",
    type == "both" & isoprene <= 30 & monoterpenes <= 5.1 ~ "medium-medium",
    type == "both" & isoprene <= 30 ~ "medium-high",
    type == "both" & monoterpenes < 2 ~ "high-low",
    type == "both" & monoterpenes <= 5.1 ~ "high-medium",
    type == "both" ~ "high-high",
    # For NE (not emitting) or others
    TRUE ~ NA_character_
  ))
  
  result <- result #|>  group_by(type, sub_type) |> dplyr::summarise(n = dplyr::n())
  return(result)
  
  
}
PieDonut_custom_colors_v2 <- function(data, 
                                      type_colors = NULL,
                                      subtype_colors = NULL,
                                      type_labels = NULL,
                                      subtype_labels = NULL,
                                      subtype_order = NULL,
                                      title = "Emission Types Distribution",
                                      show_counts = TRUE,
                                      use_repel = FALSE,
                                      inner_radius = 0.3,
                                      middle_radius = 0.7,
                                      outer_radius = 1.0,
                                      label_threshold = 0.02,
                                      segment_length = 0.15,
                                      repel_segment_length = 0.5,
                                      segment_color = "gray50",
                                      segment_size = 0.5) {
  
  # Prepare data
  df <- data %>%
    ungroup() %>%
    mutate(
      sub_type = ifelse(is.na(sub_type), "total", sub_type)
    )
  
  # Calculate totals for each type (inner ring)
  type_totals <- df %>%
    group_by(type) %>%
    summarise(total = n()) %>%
    ungroup() %>%
    mutate(
      proportion = total / sum(total),
      end_angle = cumsum(proportion) * 2 * pi,
      start_angle = lag(end_angle, default = 0),
      mid_angle = (start_angle + end_angle) / 2,
      segment_id = paste0(type, "_total")
    )
  
  # Apply type labels if provided
  if (!is.null(type_labels)) {
    type_totals <- type_totals %>%
      mutate(
        type_display = ifelse(type %in% names(type_labels), 
                              type_labels[type], 
                              type)
      )
  } else {
    type_totals <- type_totals %>%
      mutate(type_display = type)
  }
  
  # Calculate counts for subtypes
  subtype_counts <- df %>%
    group_by(type, sub_type) %>%
    summarise(count = n(), .groups = "drop")
  
  # Apply custom ordering of subtypes if specified
  if (!is.null(subtype_order)) {
    subtype_counts <- subtype_counts %>%
      group_by(type) %>%
      group_modify(~ {
        current_type <- .y$type
        if (current_type %in% names(subtype_order)) {
          .x %>% arrange(factor(sub_type, levels = subtype_order[[current_type]], ordered = TRUE))
        } else {
          .x %>% arrange(sub_type)
        }
      }) %>%
      ungroup()
  } else {
    subtype_counts <- subtype_counts %>%
      arrange(type, sub_type)
  }
  
  # Join with type angles and compute positions
  subtype_data <- subtype_counts %>%
    left_join(
      type_totals %>% select(type, type_start = start_angle, type_end = end_angle),
      by = "type"
    ) %>%
    group_by(type) %>%
    mutate(
      prop_within = count / sum(count),
      cumsum_within = cumsum(prop_within)
    ) %>%
    ungroup() %>%
    mutate(
      start_angle = type_start + (type_end - type_start) * (cumsum_within - prop_within),
      end_angle = type_start + (type_end - type_start) * cumsum_within,
      mid_angle = (start_angle + end_angle) / 2,
      segment_id = paste0(type, "_", sub_type)
    )
  
  # Apply subtype labels if provided
  if (!is.null(subtype_labels)) {
    subtype_data <- subtype_data %>%
      mutate(
        subtype_display = ifelse(sub_type %in% names(subtype_labels), 
                                 subtype_labels[sub_type], 
                                 sub_type)
      )
  } else {
    subtype_data <- subtype_data %>%
      mutate(subtype_display = sub_type)
  }
  
  # Apply type labels to subtype_data
  if (!is.null(type_labels)) {
    subtype_data <- subtype_data %>%
      mutate(
        type_display = ifelse(type %in% names(type_labels), 
                              type_labels[type], 
                              type)
      )
  } else {
    subtype_data <- subtype_data %>%
      mutate(type_display = type)
  }
  
  # Set default type colors if not provided
  if (is.null(type_colors)) {
    type_colors <- c(
      "NE" = "#999999", 
      "both" = "#E69F00", 
      "iso" = "#56B4E9", 
      "mono" = "#009E73"
    )
  }
  
  # Create color mapping for all segments
  all_segments <- unique(c(type_totals$segment_id, subtype_data$segment_id))
  color_mapping <- setNames(rep(NA, length(all_segments)), all_segments)
  
  # Assign colors to type totals (inner ring)
  for (type in names(type_colors)) {
    segment_id <- paste0(type, "_total")
    if (segment_id %in% names(color_mapping)) {
      color_mapping[segment_id] <- type_colors[type]
    }
  }
  
  # Assign colors to subtypes (outer ring)
  if (!is.null(subtype_colors)) {
    for (segment_id in names(subtype_colors)) {
      if (segment_id %in% names(color_mapping)) {
        color_mapping[segment_id] <- subtype_colors[segment_id]
      }
    }
  } else {
    for (type in unique(subtype_data$type)) {
      type_subtypes <- subtype_data %>%
        filter(type == !!type, sub_type != "total") %>%
        pull(segment_id)
      
      if (length(type_subtypes) > 0) {
        base_color <- ifelse(type %in% names(type_colors), 
                             type_colors[type], 
                             "#CCCCCC")
        base_rgb <- col2rgb(base_color) / 255
        
        for (i in seq_along(type_subtypes)) {
          factor <- 0.3 + 0.7 * (i - 1) / max(1, (length(type_subtypes) - 1))
          subtype_color <- rgb(
            base_rgb[1] * factor,
            base_rgb[2] * factor,
            base_rgb[3] * factor
          )
          color_mapping[type_subtypes[i]] <- subtype_color
        }
      }
    }
  }
  
  # Prepare outer ring data (without "total")
  outer_ring_data <- subtype_data %>%
    filter(sub_type != "total")
  
  # Build base plot
  p <- ggplot() +
    geom_arc_bar(
      data = type_totals,
      aes(
        x0 = 0, y0 = 0,
        r0 = inner_radius,
        r = middle_radius,
        start = start_angle,
        end = end_angle,
        fill = segment_id
      ),
      color = "white",
      size = 0.5
    ) +
    geom_arc_bar(
      data = outer_ring_data,
      aes(
        x0 = 0, y0 = 0,
        r0 = middle_radius,
        r = outer_radius,
        start = start_angle,
        end = end_angle,
        fill = segment_id
      ),
      color = "white",
      size = 0.5
    ) +
    theme_void() +
    theme(
      plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
      legend.position = "none"
    ) +
    labs(title = title) +
    scale_fill_manual(values = color_mapping)
  
  # Add labels if requested
  if (show_counts) {
    # Inner ring labels (with newline for readability)
    inner_labels <- type_totals %>%
      mutate(
        label = paste0(type_display, "\n(", total, ")"),
        label_x = (inner_radius + middle_radius) / 2 * sin(mid_angle),
        label_y = (inner_radius + middle_radius) / 2 * cos(mid_angle)
      )
    
    p <- p +
      geom_text(
        data = inner_labels,
        aes(x = label_x, y = label_y, label = label),
        size = 3.5,
        fontface = "bold"
      )
    
    # Prepare outer ring labels data
    outer_labels <- outer_ring_data %>%
      mutate(
        # Count on same line as subtype label
        label = paste0(subtype_display, " (", count, ")"),
        # Only show labels for segments above threshold
        show_label = count / sum(outer_ring_data$count) >= label_threshold
      ) %>%
      filter(show_label) %>%
      mutate(
        # Point on the outer edge of donut
        point_x = outer_radius * sin(mid_angle),
        point_y = outer_radius * cos(mid_angle),
        
        # Text alignment for non-repel case
        hjust = ifelse(sin(mid_angle) > 0, 0, 1),
        vjust = 0.5
      )
    
    if (use_repel) {
      # Check if ggrepel is installed
      if (!requireNamespace("ggrepel", quietly = TRUE)) {
        warning("Package 'ggrepel' is not installed. Falling back to standard labels.")
        use_repel <- FALSE
      } else {
        # For repel, place labels at a fixed distance (repel_segment_length) from outer edge
        outer_labels <- outer_labels %>%
          mutate(
            label_x_init = (outer_radius + repel_segment_length) * sin(mid_angle),
            label_y_init = (outer_radius + repel_segment_length) * cos(mid_angle)
          )
        
        p <- p +
          ggrepel::geom_text_repel(
            data = outer_labels,
            aes(x = point_x, y = point_y, label = label),
            nudge_x = outer_labels$label_x_init - outer_labels$point_x,
            nudge_y = outer_labels$label_y_init - outer_labels$point_y,
            segment.color = segment_color,
            segment.size = segment_size,
            min.segment.length = 0,  # Always draw segment
            size = 3,
            lineheight = 0.8,
            box.padding = 0.5,
            point.padding = 0.2,
            force = 2,               # Stronger force to keep near initial radius
            max.iter = 2000,
            direction = "both"
          )
      }
    }
    
    # If not using repel, add manual segments and labels with segment_length
    if (!use_repel) {
      outer_labels <- outer_labels %>%
        mutate(
          label_x_init = (outer_radius + segment_length + 0.05) * sin(mid_angle),
          label_y_init = (outer_radius + segment_length + 0.05) * cos(mid_angle)
        )
      
      p <- p +
        geom_segment(
          data = outer_labels,
          aes(x = point_x, y = point_y,
              xend = label_x_init, yend = label_y_init),
          color = segment_color,
          size = segment_size
        ) +
        geom_text(
          data = outer_labels,
          aes(x = label_x_init, y = label_y_init, 
              label = label, 
              hjust = hjust, vjust = vjust),
          size = 3,
          lineheight = 0.8
        )
    }
    
    # Compute limits to include all elements
    all_coords <- c(
      outer_labels$point_x, outer_labels$point_y,
      if (use_repel) outer_labels$label_x_init else outer_labels$label_x_init,
      inner_labels$label_x, inner_labels$label_y
    )
    max_abs <- max(abs(all_coords), na.rm = TRUE)
    plot_limit <- max_abs * 1.1
    
    p <- p + coord_fixed(
      xlim = c(-plot_limit, plot_limit),
      ylim = c(-plot_limit, plot_limit),
      clip = "off"
    )
    
  } else {
    p <- p + coord_fixed(
      xlim = c(-outer_radius, outer_radius) * 1.1,
      ylim = c(-outer_radius, outer_radius) * 1.1,
      clip = "off"
    )
  }
  
  return(p)
}





#' Title
#'
#' @param df all compounds, all samples sdt
#'
#' @returns plots for the 6 species emitters of both, might need to make it nicer
#' @export
#'
#' @examples
plot_species_both_internal_tradeoff <- function(df, save_path = here::here("figures", "screening", "internal_trade_off_both.png")) {
  
  #df = compounds_samples_spagg_to_keep
  # Load required libraries
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  
  # Species code to full name mapping
  species_names <- c(
    "ptre" = "Populus tremula",
    "rcat" = "Rhamnus cathartica",
    "sele" = "Salix eleagnos",
    "spen" = "Salix pentandra",
    "spur" = "Salix purpurea",
    "scin" = "Salix cinerea"
  )
  
  # Ensure sample columns are numeric
  df[, 5:ncol(df)] <- lapply(df[, 5:ncol(df)], as.numeric)
  
  # Reshape to long format, extract species code, and add full species name
  df_long <- df %>%
    select(class, 5:ncol(df)) %>%
    pivot_longer(-class, names_to = "sample", values_to = "value") %>%
    mutate(species_code = substr(sample, 1, 4)) %>%
    filter(species_code %in% names(species_names)) %>%
    mutate(species = species_names[species_code])
  
  # Sum by sample and class, then widen to get Isoprene and Monoterpene columns
  df_sum <- df_long %>%
    group_by(sample, species, class) %>%
    summarise(sum = sum(value, na.rm = TRUE), .groups = "drop") %>%
    pivot_wider(names_from = class, values_from = sum)
  
  # Create the plot
  p <- ggplot(df_sum, aes(x = Isoprene, y = `Oxygenated-monoterpenes` + Monoterpenes)) +
    geom_point(alpha = 0.7, size = 2, color = "black") +          # nicer points
    facet_wrap(~species, scales = "free") +                            # use full names
    labs(
      x = expression("Isoprene EF (" * mu * "g .g"^{-1} * " .h"^{-1} * ")"),
      y = expression("Monoterpenes EF (" * mu * "g .g"^{-1} * " .h"^{-1} * ")"),
      title = "Trade-off between isoprene and monoterpene emissions"
    ) +
    theme_bw(base_size = 12) +                                         # clean theme, larger text
    theme(                # subtle facet background
      strip.text = element_text(face = "italic"),                      # species names in italics
      panel.grid.minor = element_blank(),                              # remove minor grid lines
      plot.title = element_text(hjust = 0.5)                           # center title
    )
  
  # Sauvegarde si un chemin est fourni (non NULL)
  if (!is.null(save_path)) {
    # Créer le dossier si nécessaire
    dir.create(dirname(save_path), showWarnings = FALSE, recursive = TRUE)
    ggsave(filename = save_path, plot = p, width = 8, height = 6, dpi = 300)
    message("Graphique sauvegardé dans : ", save_path)
  }
  
  
  return(p)
}

