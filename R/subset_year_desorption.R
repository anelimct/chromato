#' Subset samples according to a choosen year and creates batches according to desorption date
#'
#' @param data , data frame qui contient le nom des tubes et leurs date de désorption
#' @param year , year of choice (2023, 2024, 2025) in a character "" format
#' @param alcanes , listes des fichiers alcanes disponibles
#'
#' @return
#' @export
#'
#' @examples
subset_year<- function(data, year, alcanes) {
  # Filtrer = garder seulement les lignes de l'année choisie
  data <- data |>
    dplyr::filter(stringr::str_detect(ID, paste0(year)))
  
  # Conversion de la date en format Date compréhensible pour Rstudio, jour/mois/année, puis regrouper par batch = désorbé la même semaine
  data <- data |> 
    dplyr::mutate(Date_desorption = as.Date(Date_desorption, format = "%d/%m/%Y")) |>  
    dplyr::mutate(    batch = dplyr::case_when(
      !startsWith(ID, "B_") ~ paste0("w", week(Date_desorption), "_", year(Date_desorption)),
      TRUE ~ NA_character_))
  
  
  blanks <- data |> 
    dplyr::filter(startsWith(ID, "B_")) |> 
    dplyr::select(Date, ID) 
  
  data_b <- data |> 
    dplyr::left_join(blanks, by = "Date", suffix = c("", "_blank"), relationship =
                       "many-to-many")|> 
    dplyr::filter(!startsWith(ID, "B_")) |> 
    dplyr::group_by(ID) |> dplyr::summarise(list(ID_blank))
  
  data <- data |> dplyr::left_join(data_b, by = "ID")
  
  
  
  # Conversion des dates des alcanes en format Date compréhensible pour Rstudio, jour/mois/année
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
