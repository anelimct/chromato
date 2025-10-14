var_paradise <- function(bvocs, paradise_reports_list, paradise_reports_iso_list, calib_files){
  #A Named list of samples names by batch
  samples_list <- samples_paradise(paradise_reports_list, calib_files) |> unlist() 
  samples_list_iso <- samples_paradise(paradise_reports_iso_list, calib_files) |> unlist()
  
  bvocs <- bvocs |> dplyr::mutate( paradise_mono = ifelse(ID %in% samples_list, T, F))|> dplyr::mutate( paradise_iso = ifelse(ID %in% samples_list_iso, T, F))
  
}

calibration_par <- function(data){
  cols_PAR <- c("PAR.début.1",        "PAR.début.2",       "PAR.milieu.1",     
                "PAR.milieu.2",  "PAR.milieu.3", "PAR.milieu.4",  "PAR.milieu.5", "PAR.milieu.6", "PAR.fin.1",     "PAR.fin.2") 
  ##corection des par du début a 23/05/2025 inclu 74
  
  soixante_quatorze <- interval( ymd_h("2023-06-15 00"),  ymd_h("2025-05-23 00"))
    
  ##corection des par du 03/06/2025 30.1
  
  
  trente <- interval( ymd_h("2025-06-03 00"),  ymd_h("2025-06-04 00"))

  
  ##after no correction
  
  data <- data |> 
    dplyr::mutate(across(all_of(cols_PAR), as.numeric)) |> 
    dplyr::mutate(across(all_of(cols_PAR), 
                         ~ ifelse(.data$interval %within% soixante_quatorze, 
                                  1.934168 * .x + 0.283437, 
                                  .x)))
  
  
  data <- data |> 
    dplyr::mutate(across(all_of(cols_PAR), as.numeric)) |> 
    dplyr::mutate(across(all_of(cols_PAR), 
                         ~ ifelse(.data$interval %within% trente, 
                                  4.734 * .x + 12.560, 
                                  .x)))

}

# selectionner tous les samples qui sont perdu ou erreur 

filter_out_samples <- function(data, T_max, µ_max_T, type, compound ){
  # selectionner tous les samples qui sont perdu ou erreur 
  #selectionner tous les samples pour les quelques on a perdu les ibutton values in 
  #tous les samples qui on une moyenne de température supérieur à 40 degres (faire une fonction avec un paramètre de,température à ne pas dépasser pour la moyenne ou pour les valeurs extremes)
  #Tous les échantillons qui ne ce sont pas ouverts dans paradise
  cols_PAR <- c("PAR.début.1",        "PAR.début.2",       "PAR.milieu.1",     
                "PAR.milieu.2",  "PAR.milieu.3", "PAR.milieu.4",  "PAR.milieu.5", "PAR.milieu.6", "PAR.fin.1",     "PAR.fin.2")
  blancs <- data |> dplyr::filter(Taxon == "BLANC") |>  dplyr::mutate(mean_PAR = rowMeans(dplyr::across( all_of(cols_PAR)), na.rm = TRUE)) |>  dplyr::mutate(issue = FALSE)
  data <- data |> dplyr::filter(Taxon != "BLANC") |>  dplyr::mutate(mean_PAR = rowMeans(dplyr::across( all_of(cols_PAR)), na.rm = TRUE))
  
  
  
  
  if(compound == "iso"){
    data <- data |> 
      dplyr::mutate( issue = dplyr::if_else(
        Etat %in% c("perdu", "Perdu", "Erreur") |
          purrr::map_lgl(values_T_in, ~ any(.x > T_max)) |
          purrr::map_dbl(values_T_in, ~ mean(.x)) > µ_max_T |
          purrr::map_dbl(values_T_in, ~ mean(.x)) < 20 |
          purrr::map_lgl(values_T_in, ~ length(.x) == 0) | mean_PAR < 500| #496
          paradise_iso ==FALSE, TRUE, FALSE))
    
  } else {
    data <- data |> 
      dplyr::mutate( issue = dplyr::if_else(
        Etat %in% c("perdu", "Perdu", "Erreur") |
          purrr::map_lgl(values_T_in, ~ any(.x > T_max)) |
          purrr::map_dbl(values_T_in, ~ mean(.x)) > µ_max_T |
          purrr::map_lgl(values_T_in, ~ length(.x) == 0) | mean_PAR < 500 | #496
          paradise_mono ==FALSE, TRUE, FALSE
      ))
  }
  

  
  if(type == "failed"){
    data <- data |> dplyr::group_by(ID_pairs) |> 
      dplyr::filter( all(issue)) |> 
      dplyr::distinct(ID_pairs, .keep_all = TRUE)
    
  } else {
    data <- data |>  
      dplyr::filter( ! issue) |> rbind(blancs)
  }
  
 return(data)
 
}


paired_samples <- function (bvocs){
  
  data <- bvocs |> dplyr::mutate(ID_pairs = dplyr::if_else(stringr::str_detect(ID, "2023"),  stringr::str_remove(ID, "_[^_]*_"), ID))
                                   
  
  return(data)
}




summarize_lost_samples <- function(data, T_max, µ_max_T) {
  data <- data |> dplyr::filter(Taxon != "BLANC")
  
  
  
  data <- data |> 
    dplyr::mutate( issue = dplyr::if_else(
      Etat %in% c("perdu", "Perdu", "Erreur") |
        purrr::map_lgl(values_T_in, ~ any(.x > T_max))|
        purrr::map_dbl(values_T_in, ~ mean(.x)) > µ_max_T |
        purrr::map_lgl(values_T_in, ~ length(.x) == 0) |
        paradise==FALSE, TRUE, FALSE
    ))
  
  total_counts <- data |>
    dplyr::summarize(total_count = dplyr::n(), .groups = 'drop')
  
  # Create a summary table of lost samples
  lost_samples_summary <- data |>
    dplyr::filter(issue) |>
    dplyr::group_by(ID) |> 
    dplyr::summarize(
      lost_count = dplyr::n(),
      lost_reasons = toString(unique(
        dplyr::case_when(
          Etat %in% c("perdu", "Perdu", "Erreur") ~ "Status: Perdu/Erreur",
          purrr::map_lgl(values_T_in, ~ any(.x > T_max)) ~ "Temp: Exceeds Max",
          purrr::map_dbl(values_T_in, ~ mean(.x)) > µ_max_T ~ "Avg Temp: Exceeds Max",
          purrr::map_lgl(values_T_in, ~ length(.x) == 0) ~ "Temp Values: Missing",
          paradise == FALSE ~ "Paradise: Not Opened",
          TRUE ~ ""
        )
      ))
    )
  # Join the total counts and lost counts
  summary_table <- lost_samples_summary |>
    dplyr::mutate(lost_count = ifelse(is.na(lost_count), 0, lost_count),
           lost_reasons = ifelse(is.na(lost_reasons), "None", lost_reasons)) |> 
    dplyr::filter(lost_count != 0)
  
  return(summary_table)
}

summarize_field <- function(data_iso, data_mono) {
  data_iso <-  data_iso |> dplyr::filter(Taxon != "BLANC") |> 
    dplyr::group_by(gragg, Taxon) |>
    dplyr::summarise(
      distinct_origins_field_iso = dplyr::n_distinct(paste(stringr::str_sub(Latitude..WGS84., 1, 4), stringr::str_sub(Longitude..WGS84., 1, 4), sep = "_")), 
      n_entries_iso = dplyr::n()
      ) |>
    dplyr::ungroup()
  
  data_mono <- data_mono |> dplyr::filter(Taxon != "BLANC") |> 
    dplyr::group_by(gragg, Taxon) |>
    dplyr::summarise(
      distinct_origins_field_mono = dplyr::n_distinct(paste(stringr::str_sub(Latitude..WGS84., 1, 4), stringr::str_sub(Longitude..WGS84., 1, 4), sep = "_")),
      n_entries_mono = dplyr::n()) |>
    dplyr::ungroup()
  
  data <-  dplyr::left_join(data_mono, data_iso)
}


  
