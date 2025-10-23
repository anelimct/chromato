

fusionner_listes <- function(ma_liste) {
  # Sélectionner les colonnes et fusionner en une étape
  ma_liste |> 
    # Sélectionner les colonnes désirées
    purrr::map(~ {
      colonnes_samples <- names(.x)[grepl("^[a-zA-Z]{4}_", names(.x))]
      
      .x |> 
        dplyr::select(
          dplyr::any_of( c("compound", "smiles", "class", "inchikey")),
          dplyr::all_of(colonnes_samples)
        )
    }) |> 
    # Fusionner tous les dataframes
   purrr::reduce(dplyr::full_join, by = c("compound", "smiles", "class", "inchikey"))
    }





analyser_statut_compounds <- function(df) {
  # Identifier les colonnes samples (qui commencent par 4 lettres + _)
  colonnes_samples <- names(df)[grepl("^[a-zA-Z]{4}_", names(df))]
  
  df <- df |>   dplyr::mutate(dplyr::across(tidyselect::everything(), as.character))
  
  # Créer un résumé par compound
  resume_compounds <- df |> 
    tidyr::pivot_longer(
      cols = all_of(colonnes_samples),
      names_to = "sample",
      values_to = "valeur"
    ) |> 

    dplyr::group_by(compound) |> 
    dplyr::summarise(
      # Statuts de détection
      nb_total = n(),
      nb_nd = sum(valeur == "nd", na.rm = TRUE),
      nb_tr = sum(valeur == "tr", na.rm = TRUE),
      nb_quantifies = sum(!valeur %in% c("nd", "tr") & !is.na(valeur)), 
      categorie = case_when(
        nb_quantifies > 0 ~ "Quantifié",
        nb_tr > 0 ~ "Traces uniquement",
        TRUE ~ "Non détecté"
      ),
      
      .groups = "drop"
    )
}





trouver_top_samples <- function(table_statut, table_quantites) {
  
  compounds_cibles <- table_statut |> 
    dplyr::filter(categorie %in% c("Traces uniquement", "Non détecté")) %>%
    dplyr::pull(compound)
  
  # Préparer la table de quantités (garder seulement les compounds cibles)
  quantites_filtrees <- table_quantites |> 
    dplyr::filter(compound %in% compounds_cibles) 
    
    colonnes_samples <- names(quantites_filtrees)[grepl("^[a-zA-Z]{4}_", names(quantites_filtrees))]
  
  
    quantites_filtrees <-quantites_filtrees |> 
    tidyr::pivot_longer(
      cols = tidyselect::all_of(colonnes_samples),
      names_to = "sample",
      values_to = "quantity_ng"
    ) |> 
    # Grouper par compound et trouver les 2 plus grandes aires
    dplyr::group_by(compound) |> 
    dplyr::arrange(desc(quantity_ng), .by_group = TRUE)  |> 
    slice_head(n = 2) |> 
    # Créer une colonne pour le rang
    dplyr::mutate(rang = paste0("top_", row_number())) |>
    ungroup() |> 
    # Remettre en format wide
    dplyr::select(compound, rang, sample, quantity_ng) |> 
    tidyr::pivot_wider(
      names_from = rang,
      values_from = c(sample,quantity_ng),
      names_glue = "{rang}_{.value}"
    )
  
  # Fusionner avec les informations de statut
  resultat_final <- table_statut  |> 
    dplyr::filter(categorie %in% c("Traces uniquement", "Non détecté"))  |> 
    dplyr::left_join(quantites_filtrees, by = "compound")
  
  return(resultat_final)
}





merge_datasets <- function(compounds_table, bvocs_samples, valid_samples_mono, paradise_reports_mono_ER){
  #samples pour lesquels les mono ont été analysés (.cdf)
  rows_to_plot <- names( fusionner_listes(paradise_reports_mono_ER))[grepl("^[a-zA-Z]{4}_", names( fusionner_listes(paradise_reports_mono_ER)))] 
  #samples qui repectect les condition <43 °C , PAR min (pas.cdf avec des majuscules)
  rows_to_keep <- valid_samples_mono$ID
  
  df <- compounds_table
  colonnes_samples <- names( df)[grepl("^[a-zA-Z]{4}_", names( df))]
  
  df <- compounds_table |> dplyr::group_by(class) |> dplyr::summarise(dplyr::across(tidyselect::all_of(colonnes_samples), ~ sum(as.numeric(.x), na.rm = TRUE),
             .names = "{.col}")) |> dplyr::ungroup() |>    tidyr::pivot_longer(cols = -class,names_to = "sample",values_to = "value"
    ) |> tidyr::pivot_wider(names_from = class,values_from = value) |> 
    dplyr::filter(sample %in% rows_to_plot) |> dplyr::mutate(ID = stringr::str_replace(sample, ".cdf", "")) |> dplyr::mutate(ID = stringr::str_to_upper(ID)) 
    
  cols_PAR <- c("PAR.début.1",        "PAR.début.2",       "PAR.milieu.1",     
                "PAR.milieu.2",  "PAR.milieu.3", "PAR.milieu.4",  "PAR.milieu.5", "PAR.milieu.6", "PAR.fin.1",     "PAR.fin.2")
  
  bvocs_samples_ER <- merge(bvocs_samples, df, by = "ID") |> dplyr::filter(!is.na(Leaves_DM)) |>    dplyr::mutate(across(all_of(cols_PAR), as.numeric)) |>  dplyr::mutate(PAR_algo = rowMeans(dplyr::across( all_of(cols_PAR)), na.rm = TRUE)) |> 
    dplyr::mutate(mean_T = purrr::map_dbl(values_T_in, ~ mean(.x))) |> dplyr::mutate(T_algo_K = mean_T + 273,15) |> 
    dplyr::select(ID, Monoterpenes, Isoprene, Taxon, PAR_algo, T_algo_K ) |>  tidyr::pivot_longer(cols = c(Isoprene, Monoterpenes),names_to = "Compound",values_to = "Emission") |> dplyr::mutate(Standardized = FALSE) |> dplyr::mutate(Compound = stringr::str_to_lower(Compound))
  
  data_mono <- bvocs_samples_ER |>  dplyr::filter(Compound == "Monoterpenes")
  data_iso <- bvocs_samples_ER |>  dplyr::filter(Compound == "Isoprene")
  
  plot_mono <- ggplot(data_mono) + 
    aes(x = T_algo_K - 273,15, y = Emission)+
    geom_point(color = "blue", size = 2, alpha = 0.7)
  
  plot_iso <- ggplot(data_iso) + 
    aes(x = T_algo_K - 273,15, y = Emission)+
    geom_point(color = "blue", size = 2, alpha = 0.7)
  
  bvocs_samples_EF <-bvocs_samples_ER |> dplyr::filter( ID %in% rows_to_keep ) |>  apply_standardization() |>  dplyr::mutate(EF = ES_iso_G93) |>  dplyr::select(EF, Taxon, Compound) |> dplyr::mutate(source = "field") |>  dplyr::mutate(Taxon = stringr::str_replace_all(Taxon, " ", "_")) |> dplyr::mutate(EF = EF/1000)
  
  
}



      
      
      