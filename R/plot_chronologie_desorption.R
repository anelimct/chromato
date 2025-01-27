
#' Creates a visual mapping of desorption date in the year_subplot
#'
#' @param data  a dataframe qui contient le nom des tubes fait la même année et l'attribution a des batches
#' @return a plot with time in the x axis
#' @export
#'
#' @examples
chronologie <- function(data) {
  summary_data <- data |>
    dplyr::filter(!stringr::str_detect(batch, "NA")) |>
    dplyr::group_by(batch) |>
    dplyr::summarise(
      first_date = min(Date_desorption, na.rm = TRUE),
      last_date = max(Date_desorption, na.rm = TRUE),
      num_samples = dplyr::n(),
      num_exploitables = sum(Pics_exploitables != "Non", na.rm = TRUE),  # Compte des exploitables
      num_paradise = sum(Remarque_pics == "OK", na.rm = FALSE),  # Compte des exploitables paradise
      .groups = "drop"
    )
  
  # Préparer les données pour inclure les trois types
  bar_data <- summary_data |>
    tidyr::pivot_longer(
      cols = c(num_samples, num_exploitables, num_paradise),
      names_to = "type",
      values_to = "count"
    ) |>
    dplyr::mutate(
      y_position = dplyr::case_when(
        type == "num_samples" ~ as.numeric(factor(batch)),  # Position centrale (barre rose)
        type == "num_exploitables" ~ as.numeric(factor(batch)) - 0.2,  # Barre verte, décalée vers le bas
        type == "num_paradise" ~ as.numeric(factor(batch)) -0.4   # Barre bleue, décalée vers le haut
      ),
      color = dplyr::case_when(
        type == "num_samples" ~ "Total",
        type == "num_exploitables" ~ "Exploitables",
        type == "num_paradise" ~ "Paradise"
      )
    )
  
  # Créer le graphique
  p <- ggplot() +
    # Barres avec légende
    geom_segment(data = bar_data, 
                 aes(x = first_date, xend = last_date, 
                     y = y_position, yend = y_position, color = color), 
                 linewidth = 4, alpha = 0.7) +
    # Annotations
    geom_text(data = bar_data,
              aes(x = last_date + 2, y = y_position, label = count, color = color),
              hjust = 0, size = 3.5) +
    labs(
      title = "Événements de désorption par batch",
      x = "Date",
      y = "Batch",
      color = "Type d'échantillon",  # Titre de la légende
      caption = "La longueur de la barre représente une durée, le nombre d'échantillons est indiqué en chiffre"
    ) +
    scale_y_continuous(
      breaks = 1:length(unique(summary_data$batch)), 
      labels = unique(summary_data$batch)
    ) +
    scale_color_manual(
      values = c("Total" = "#ab4974", "Exploitables" = "#4caf50", "Paradise" = "#85d3e2")
    ) +
    theme_minimal() +
    theme(
      axis.text.y = element_text(size = 10),
      axis.title.y = element_text(size = 12),
      plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
      legend.title = element_text(size = 12),
      legend.text = element_text(size = 10)
    )
  
  return(p)
}



