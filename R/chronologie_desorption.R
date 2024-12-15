
#' Creates a visual mapping of desorption date in the year_subplot
#'
#' @param data 
#'
#' @return a plot with time in the x axis
#' @export
#'
#' @examples
chronologie <- function(data){
  summary_data <- data  |>
    dplyr::filter(!stringr::str_detect(batch, "NA")) |> 
    dplyr::group_by(batch) |> 
    dplyr::summarise(
    first_date = min(Date_desorption, na.rm = TRUE),
    last_date = max(Date_desorption, na.rm = TRUE),
    num_samples = dplyr::n(),
    .groups = "drop"
  )
  
  # # Créer le graphique
   p <- ggplot(summary_data, aes(x = first_date, xend = last_date, y = factor(batch), yend = factor(batch))) +
     geom_segment(linewidth = 2, color = "#ab4974") +  # Ligne de chaque batch
     geom_text(aes(x = first_date, label = num_samples), hjust = -0.5, vjust = -1, color = "black") +  # Nombre d'échantillons
     labs(
      title = "Événements de désorption par batch",
       x = "Date",
       y = "Batch",
       caption = "Ligne = durée de désorption, Nombre = échantillons par batch"
     ) +
    theme_minimal() +
     theme(axis.text.y = element_text(size = 10),
           axis.title.y = element_text(size = 12))
  return(p)
}


