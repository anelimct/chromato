#' Subset samples according to a choosen year and creates batches according to desorption date
#'
#' @param data , data frame
#' @param year , year of choice (2023, 2024, 2025) in a character "" format
#' @param alcanes 
#'
#' @return
#' @export
#'
#' @examples
subset_year<- function(data, year, alcanes) {
  # Filtrer pour les lignes de 2024
  data <- data |>
    dplyr::filter(stringr::str_detect(ID, paste0(year)))
  
  # Conversion en format Date
  data <- data |> 
    dplyr::mutate(Date_desorption = as.Date(Date_desorption, format = "%d.%m.%Y")) |>  
    dplyr::mutate(batch = paste0( "w", lubridate::week(Date_desorption), "_", year))
  
  # Conversion des dates des alcanes en format Date
  alcanes <- alcanes |> 
    dplyr::mutate(date = as.Date(date, format = "%d/%m/%Y"))

  # Associer à chaque ligne de 'data' le fichier alcane le plus proche
  data <- data |>
    dplyr::rowwise() |>
    dplyr::mutate(
      alcane_closest = list(alcanes$new_name[which.min(abs(difftime(Date_desorption, alcanes$date, units = "days")))])
    ) |>
    dplyr::ungroup()
  
  return(data) 
}
