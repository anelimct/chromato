#' Subset samples according to a choosen year and creates batches according to desorption date
#'
#' @param data , data frame
#' @param year , year of choice (2023, 2024, 2025) in a character "" format
#'
#' @return
#' @export
#'
#' @examples
subset_year<- function(data, year) {
  # Filtrer pour les lignes de 2024
  data <- data |>
    dplyr::filter(stringr::str_detect(ID, paste0(year)))
  
  # Conversion en format Date
  data <- data |> 
    dplyr::mutate(Date_desorption = as.Date(Date_desorption, format = "%d.%m.%Y")) |>  
    dplyr::mutate(batch = paste0( "w", lubridate::week(Date_desorption), "_", year))
  
  return(data)
}
