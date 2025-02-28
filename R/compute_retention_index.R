imput_missing_alcanes <- function(alcanes) {
  # Calculer la différence moyenne entre "undecane" et "decane"
  mean_c10_c11 <- mean(alcanes$undecane - alcanes$decane, na.rm = TRUE)
  
  # Imputer les valeurs manquantes pour "undecane" et "decane"
  alcanes$undecane <- ifelse(is.na(alcanes$undecane) & !is.na(alcanes$decane),
                             alcanes$decane + mean_c10_c11,
                             alcanes$undecane)
  alcanes$decane <- ifelse(is.na(alcanes$decane) & !is.na(alcanes$undecane),
                           alcanes$undecane - mean_c10_c11,
                           alcanes$decane)
  
  # Calculer la différence moyenne entre "pentadecane" et "hexadecane"
  mean_c15_c16 <- mean(alcanes$hexadecane - alcanes$pentadecane, na.rm = TRUE)
  
  # Imputer les valeurs manquantes pour "hexadecane" et "pentadecane"
  alcanes$hexadecane <- ifelse(is.na(alcanes$hexadecane) & !is.na(alcanes$pentadecane),
                               alcanes$pentadecane + mean_c15_c16,
                               alcanes$hexadecane)
  alcanes$pentadecane <- ifelse(is.na(alcanes$pentadecane) & !is.na(alcanes$hexadecane),
                                alcanes$hexadecane - mean_c15_c16,
                                alcanes$pentadecane)
  return(alcanes)
}



calculate_IR_exp_for_one_column <- function(samples_RT_alc, paradise_report, sample_col) {
  paradise_report <- paradise_report |>
    dplyr::rowwise() |>
    dplyr::mutate(
      !!sample_col := {
        ligne <- samples_RT_alc |> dplyr::filter(ID == sample_col)
        
        if (nrow(ligne) == 0) {
          stop(paste("No matching ID found in samples_RT_alc for:", sample_col))
        } 
          
        RT_n <- dplyr::case_when(
          nb_C_n == 10 ~ ligne$decane,
          nb_C_n == 15 ~ ligne$pentadecane
        )
        RT_n_s <- dplyr::case_when(
          nb_C_n_s == 11 ~ ligne$undecane,
          nb_C_n_s == 16 ~ ligne$hexadecane
        )
        100 * (nb_C_n + (`Est. Retention Time (min)` - RT_n) / (RT_n_s - RT_n))
      }
    )
  return(paradise_report)
}

calculate_IR_exp_for_all_columns <- function(samples_RT_alc, paradise_report, sample_cols) {
  for (i  in 1:length(sample_cols)) {
    sample_col <- sample_cols[i]
    print(sample_col)
    # Call the function for the current column
    paradise_report <- calculate_IR_exp_for_one_column(samples_RT_alc, paradise_report, sample_col)
  }
  
  
  
  combined_df <- paradise_report |>
    dplyr::rowwise() |>
    dplyr::mutate(RI_exp = list(unique(dplyr::c_across(all_of(sample_cols))))) |>
    dplyr::mutate(range_RI = list(sapply(RI_exp, function(x) abs(x - RI))))|>
    dplyr::ungroup() |>
    dplyr::select("Compound Name", "Match Quality", "Est. Retention Time (min)", "Hit 1: Probability", "Hit 1: CAS", "New_CAS_Name", "RI", "nb_C_n", "nb_C_n_s", "RI_exp", "range_RI")
  
  return(combined_df)
}



process_paradise_reports <- function(paradise_reports_list, samples_list, samples_RT_alc) {
  # Itérer sur chaque élément des listes nommées
  purrr::map2(paradise_reports_list, samples_list, ~ {
    paradise_report <- .x
    sample_cols <- .y
    updated_report <- calculate_IR_exp_for_all_columns(samples_RT_alc, paradise_report, sample_cols)
    return(updated_report)
  })
}




compute_retention_index <- function(alcanes, bvocs, paradise_reports_list, calib_files) {
  
  #dataframe with samples ID as rows and columns undecane, decane, pentane et hexadecane RT
  samples_RT_alc <- dplyr::left_join(bvocs, alcanes, by = c("alcane_closest" = "new_name")) |> dplyr::mutate( ID = paste0(ID, ".CDF")) 
  
  #A Named list of samples names by batch
  samples_list <- samples_paradise(paradise_reports_list, calib_files) |> purrr::map( ~ paste0(.x, ".CDF"))
  
  RI <- process_paradise_reports(paradise_reports_list, samples_list, samples_RT_alc)
  
  return(RI)  

}






