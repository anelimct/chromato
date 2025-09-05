

get_compound_stats <- function(paradise_reports_mono_ER_list) {
  # 1. Combiner tous les rapports en un seul dataframe
  combined_data <- dplyr::bind_rows(paradise_reports_mono_ER_list) |> 
    dplyr::select(New_CAS_Name, smiles, matches("\\.CDF$")) |>  # Sélectionne les colonnes .CDF
    dplyr::select(-starts_with("B_"), -contains("calib")) |>    # Exclut les blanks et calibrations
    tidyr::pivot_longer(
      cols = -c(New_CAS_Name, smiles, inchikey),
      names_to = "Sample",
      values_to = "Value"
    ) |>
    dplyr::filter(!is.na(New_CAS_Name))  # Supprimer les lignes sans identifiant CAS
  
  # 2. Calculer les statistiques pour chaque composé
  compound_stats <- combined_data |>
    dplyr::group_by(New_CAS_Name, smiles) |>
    dplyr::summarise(
      total_samples = dplyr::n(),
      nd_count = sum(Value == "nd", na.rm = TRUE),
      tr_count = sum(Value == "tr", na.rm = TRUE),
      quant_count = sum(!(Value %in% c("nd", "tr")) & !is.na(Value), na.rm = TRUE),
      .groups = 'drop'
    ) |>
    dplyr::mutate(
      detection_category = dplyr::case_when(
        quant_count >= 1 ~ "quantifier",
        quant_count == 0 & tr_count >= 1 ~ "Traces",
        quant_count == 0 & tr_count == 0 ~ "Below detection limit"
        
      )
    ) |>
    dplyr::arrange(desc(quant_count), New_CAS_Name) |>
    dplyr::select(
      New_CAS_Name, smiles,
      total_samples, 
      quantified = quant_count,
      tr_count,
      nd_count,
      detection_category
    )
  
  return(compound_stats)
}
