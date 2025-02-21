var_paradise <- function(bvocs, paradise_reports_list, calib_files){
  #A Named list of samples names by batch
  samples_list <- samples_paradise(paradise_reports_list, calib_files) |> unlist() 
  
  bvocs <- bvocs |> dplyr::mutate( paradise = ifelse(ID %in% samples_list, T, F))
  
}

# selectionner tous les samples qui sont perdu ou erreur 

filter_out_samples <- function(data, T_max, µ_max_T){
  # selectionner tous les samples qui sont perdu ou erreur 
  #selectionner tous les samples pour les quelques on a perdu les ibutton values in 
  #tous les samples qui on une moyenne de température supérieur à 40 degres (faire une fonction avec un paramètre de,température à ne pas dépasser pour la moyenne ou pour les valeurs extremes)
  data <- data |> dplyr::filter(
    Etat %in% c("perdu", "Perdu", "Erreur") |
                                  purrr::map_lgl(values_T_in, ~ any(.x > T_max))|
                                  purrr::map_dbl(values_T_in, ~ mean(.x)) > µ_max_T |
                                  purrr::map_lgl(values_T_in, ~ length(.x) == 0)
                                )
  
  data <- data |> dplyr::filter(Taxon != "BLANC")
  
}


summary_filter_out_samples <- function(data, T_max, µ_max_T){

    # Ajouter des colonnes indicatrices pour chaque condition
    data <- data |>
      dplyr::mutate(
        is_perdu_erreur = Etat %in% c("perdu", "Perdu", "Erreur"),
        any_above_T_max = purrr::map_lgl(values_T_in, ~ any(.x > T_max)),
        mean_above_µ_max_T = purrr::map_dbl(values_T_in, ~ mean(.x)) > µ_max_T,
        is_empty = purrr::map_lgl(values_T_in, ~ length(.x) == 0)
      )
    
    # Filtrer les lignes qui correspondent à au moins une condition
    filtered_data <- data |>
      dplyr::filter(is_perdu_erreur | any_above_T_max | mean_above_µ_max_T | is_empty) |> 
      dplyr::filter(Taxon != "BLANC")
    
    # Compter le nombre de TRUE pour chaque condition
    summary_table <- filtered_data |>
      dplyr::summarise(
        n_perdu_erreur = sum(is_perdu_erreur),
        n_any_above_T_max = sum(any_above_T_max),
        n_mean_above_µ_max_T = sum(mean_above_µ_max_T, na.rm = TRUE),
        n_inutton_lost = sum(is_empty)
      )
    
}
  
  
