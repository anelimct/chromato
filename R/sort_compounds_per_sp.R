apply_unknown_compounds <- function(compounds_table){
  compounds_table |>
    dplyr::mutate(
      compound = ifelse(compound == "1,3,8-p-Menthatriene", 
                        "unknown monoterpene 1", compound),
      smiles = ifelse(compound == "unknown monoterpene 1", NA, smiles),
      inchikey = ifelse(compound == "unknown monoterpene 1", NA, inchikey)
    ) |> 
    dplyr::mutate(
      compound = ifelse(compound == "α-Patchoulene", 
                        "unknown sesquiterpene 1", compound),
      smiles = ifelse(compound == "unknown sesquiterpene 1", NA, smiles),
      inchikey = ifelse(compound == "unknown sesquiterpene 1", NA, inchikey)
    ) |>  
    dplyr::mutate(
      compound = ifelse(compound == "α-Dehydro-ar-himachalene", 
                        "unknown sesquiterpene 2", compound),
      smiles = ifelse(compound == "unknown sesquiterpene 2", NA, smiles),
      inchikey = ifelse(compound == "unknown sesquiterpene 2", NA, inchikey)
    )
  
  
}



standardisation_whole_compounds_table <- function(bvocs_samples, compounds_table){
  
  # Créer un format long avec tous les compounds individuels
  colonnes_samples <- names(compounds_table)[grepl("^[a-zA-Z]{4}_", names(compounds_table))]
  
  df_long <- compounds_table |> 
    tidyr::pivot_longer(
      cols = tidyselect::all_of(colonnes_samples),
      names_to = "sample",
      values_to = "Emission"
    ) |> 
    dplyr::mutate(Emission = as.numeric(Emission)) |> 
    dplyr::mutate(ID = stringr::str_replace(sample, ".cdf", "")) |> 
    dplyr::mutate(ID = stringr::str_to_upper(ID)) 
  
  cols_PAR <- c("PAR.début.1", "PAR.début.2", "PAR.milieu.1",     
                "PAR.milieu.2", "PAR.milieu.3", "PAR.milieu.4", "PAR.milieu.5", "PAR.milieu.6", "PAR.fin.1", "PAR.fin.2")
  
  
  bvocs_samples_long <- merge(bvocs_samples, df_long, by = "ID") |> 
    dplyr::filter(!is.na(Leaves_DM)) |>    
    dplyr::mutate(across(all_of(cols_PAR), as.numeric)) |>  
    dplyr::mutate(PAR_algo = rowMeans(dplyr::across(all_of(cols_PAR)), na.rm = TRUE)) |> 
    dplyr::mutate(mean_T = purrr::map_dbl(values_T_in, ~ mean(.x))) |> 
    dplyr::mutate(T_algo_K = mean_T + 273.15) |> 
    dplyr::select( compound, class, smiles, inchikey, sample,  Emission,  PAR_algo, T_algo_K) |>  
    dplyr::mutate(Standardized = FALSE)
  
  
  bvocs_samples_standardized <- bvocs_samples_long |> 
    apply_standardization() |> 
    dplyr::mutate(
      EF_ng = dplyr::case_when(
        class %in% c("Sesquiterpenes", "Oxygenated-sesquiterpenes") ~ ES_sesqui_G93,
        TRUE ~ ES_iso_G93
      )
    ) |>  dplyr::mutate(EF_µg = EF_ng / 1000) |> 
    dplyr::mutate(
      EF_ng_T = dplyr::case_when(
        class %in% c("Monoterpenes", "Oxygenated-monoterpenes") ~ ES_mono_G93,
        TRUE ~ NA
      )
    ) |>  dplyr::mutate(EF_µg_T = EF_ng_T / 1000) |> 
    dplyr::mutate(
      EF_ng_T_Bourtsou = dplyr::case_when(
        class %in% c("Monoterpenes", "Oxygenated-monoterpenes") ~ ES_mono_G93_bourtsoukidis,
        TRUE ~ NA
      )
    ) |>  dplyr::mutate(EF_µg_T_Bourtsou = EF_ng_T_Bourtsou / 1000) 
    
  
  ##Ici on garde le choix de standardisation definitif
  
  compounds_table_standardized <- bvocs_samples_standardized |> 
    dplyr::select(compound, class,smiles, inchikey, sample, EF_µg) |> 
    dplyr::group_by(compound, class, sample, smiles, inchikey) |>  
    dplyr::summarise(EF_µg = dplyr::first(EF_µg), .groups = "drop") |>
    tidyr::pivot_wider(
      names_from = sample,
      values_from = EF_µg,
    )|>
    dplyr::mutate(
      dplyr::across(
        dplyr::matches("2023|2024"),
        ~ ifelse(class %in% c("Sesquiterpenes", "Oxygenated-sesquiterpenes") & !is.na(.x), 
                 "present", 
                 .x)
      ))
  ##  ici les nd se sont transformé en NA on a perdu cette information 
  
  compounds_table_standardized_mono_T <- bvocs_samples_standardized |> 
    dplyr::select(compound, class,smiles, inchikey, sample, EF_µg_T) |> 
    dplyr::filter(class %in% c("Monoterpenes", "Oxygenated-monoterpenes")) |> 
    dplyr::group_by(compound, class, sample, smiles, inchikey) |>  
    dplyr::summarise(EF_µg_T = dplyr::first(EF_µg_T), .groups = "drop") |>
    tidyr::pivot_wider(
      names_from = sample,
      values_from = EF_µg_T,
    )
  
  #sauvegarder dans un fichier excel les tot sdt avec T pour monoterpenes, mais ici on atoujours les singletons ! il faudra les passé dans une autre fonction
  
  compounds_table_standardized_mono_T_bourtsou <- bvocs_samples_standardized |> 
    dplyr::select(compound, class,smiles, inchikey, sample, EF_µg_T_Bourtsou) |>
    dplyr::filter(class %in% c("Monoterpenes", "Oxygenated-monoterpenes")) |>
    dplyr::group_by(compound, class, sample, smiles, inchikey) |>  
    dplyr::summarise(EF_µg_T_Bourtsou = dplyr::first(EF_µg_T_Bourtsou), .groups = "drop") |>
    tidyr::pivot_wider(
      names_from = sample,
      values_from = EF_µg_T_Bourtsou,
    )
  
  
  compounds_table_standardized_mono_LetT <- bvocs_samples_standardized |> 
    dplyr::select(compound, class,smiles, inchikey, sample, EF_µg) |> 
    dplyr::filter(class %in% c("Monoterpenes", "Oxygenated-monoterpenes")) |>
    dplyr::group_by(compound, class, sample, smiles, inchikey) |>  
    dplyr::summarise(EF_µg = dplyr::first(EF_µg), .groups = "drop") |>
    tidyr::pivot_wider(
      names_from = sample,
      values_from = EF_µg,
    )
  
  
  writexl::write_xlsx( compounds_table_standardized_mono_T, path = here::here("outputs", "comparison_mono_std", "ER_mono_TG93.xlsx"))
  writexl::write_xlsx( compounds_table_standardized_mono_LetT, path = here::here("outputs", "comparison_mono_std", "ER_mono_isoG93.xlsx"))
  writexl::write_xlsx( compounds_table_standardized_mono_T_bourtsou, path = here::here("outputs", "comparison_mono_std", "ER_mono_T_Bourtsoukidis.xlsx"))
  
  return(compounds_table_standardized)
}




times_compound_per_species <- function(compounds_table, valid_samples_mono){
  cols_to_keep <- valid_samples_mono$ID |>  stringr::str_to_lower() |> stringr::str_c( ".cdf", sep="") |>
    (\(x) x[grepl("^[a-zA-Z]{4}_", x)])()
  
  df <- compounds_table |>  dplyr::select( any_of( c("compound", "smiles",  "class" , "inchikey", paste0(cols_to_keep))))

  
  # 1. Préparer les métadonnées des échantillons
  sample_metadata <- data.frame(
    sample_id = colnames(df)[5:ncol(df)],  # À partir de la 5ème colonne
    spagg = stringr::str_extract(colnames(df)[5:ncol(df)], "^[a-z]{4}") |>  stringr::str_to_upper()  # Extraire les 4 premières lettres
  )
  
  # 2. Compter le nombre total d'échantillons par espèce
  samples_per_species <- sample_metadata |> 
    dplyr::count(spagg, name = "total_samples_species")
  
  # 3. Transformer les données en format long et calculer les statistiques
  detection_stats <- df |> 
    # Sélectionner les colonnes des composés et des échantillons
    dplyr::select(compound, all_of(sample_metadata$sample_id))|>
    # Transformer en format long
    tidyr::pivot_longer(
      cols = -compound,
      names_to = "sample_id",
      values_to = "intensity"
    ) |> 
    # Ajouter les métadonnées des espèces
    dplyr::left_join(sample_metadata, by = "sample_id")|>
    # Définir ce qu'est une détection (à adapter selon vos données)
    dplyr::mutate(
      detected = !is.na(intensity) & intensity > 0 & intensity != "nd" # ici modif pour les present sesqui que ça compte bien comme oui on l'a trouvé
    )|>
    # Grouper par composé et espèce
    dplyr::group_by(compound, spagg)|>  
    dplyr::summarise(
      samples_detected = sum(detected, na.rm = TRUE),
      total_samples_in_group = dplyr::n(),
      .groups = "drop"
    )|>
    # Ajouter le nombre total d'échantillons par espèce
    dplyr::left_join(samples_per_species, by = "spagg")|>
    # Calculer le taux de détection
    dplyr::mutate(
      detection_rate = samples_detected / total_samples_species
    )|>
    # Réorganiser les colonnes pour plus de clarté
    dplyr::select(
      compound, spagg,
      samples_detected, 
      total_samples_species,
      detection_rate
    ) |>  dplyr::filter(total_samples_species >= 3)|>  dplyr::filter(detection_rate >= 0.5) |>  dplyr::mutate(spagg = stringr::str_to_lower(spagg))
}








compound_mean_sp <- function(compounds_table, valid_samples_mono, times_compound_sp){
  
  cols_to_keep <- valid_samples_mono$ID |>  stringr::str_to_lower() |> stringr::str_c( ".cdf", sep="") |>
    (\(x) x[grepl("^[a-zA-Z]{4}_", x)])()
  
  df <- compounds_table |>  dplyr::select( any_of( c("compound", "smiles",  "class" , "inchikey", paste0(cols_to_keep))))
  
  
 spagg_to_keep <- unique(times_compound_sp$spagg)
 
 ## pour toutes les colonnes qui commence par le même spagg = les memes premiere lettre et don ces premières lettres sont dans spagg to keep faire 
 ##la moyenne et l'erreur standard entre parenthèse. Pour connaitre pas combien diviser pour faire la moyenne il y a des subtilités. Déjà pour faire la moyenne d'un composé pour cet spagg il faut que ce couple spagg/ compound soit dans times_compound_sp  ou il y a une colonne spagg et une colonne compounds. 
 ## pour savoir par combien diviser, pour les class %in% ("Sesquiterpenes", "Oxygenated-sesquiterpenes") prendre seulement le nombre de samples fait en 2025
 ## pour le reste des class il faut diviser par le total_samples_species de times_compound_sp pour ce couple espèce composé 
  
 # Fonction pour extraire le préfixe spagg d'un nom de colonne
 extract_spagg <- function(col_name) {
   # Supprimer l'extension .cdf et prendre les 4 premières lettres
   base_name <- sub("\\.cdf$", "", col_name)
   substr(base_name, 1, 4)
 }
 
 # Fonction pour calculer la moyenne et l'erreur standard pour un spagg donné
 calculate_mean_se <- function(compound_row, spagg, compound_name) {
   # Trouver les colonnes correspondant au spagg
   spagg_cols <- names(compound_row)[sapply(names(compound_row), function(x) extract_spagg(x) == spagg)]
   
   if (length(spagg_cols) == 0) {
     return(NA)
   }
   
   # Extraire les valeurs pour ces colonnes
   values <- as.numeric(compound_row[spagg_cols])
   
   # Obtenir la classe du composé
   compound_class <- compound_row$class
   
   # Trouver le nombre d'échantillons pour la division
   if (compound_class %in% c("Sesquiterpenes", "Oxygenated-sesquiterpenes")) {
     # Pour ces classes, utiliser seulement les échantillons de 2025
     # Vous devrez peut-être adapter cette logique selon vos données
     n_samples <- length(spagg_cols)  # À adapter selon vos critères 2025
   } else {
     # Pour les autres classes, utiliser total_samples_species de times_compound_sp
     species_compound_info <- times_compound_sp[
       times_compound_sp$spagg == spagg & times_compound_sp$compound == compound_name,
     ]
     
     if (nrow(species_compound_info) > 0) {
       n_samples <- species_compound_info$total_samples_species[1]
     } else {
       n_samples <- length(spagg_cols)  # Valeur par défaut
     }
   }
   
   # Calculer la moyenne
   mean_val <- mean(values, na.rm = TRUE)
   
   # Calculer l'erreur standard
   se_val <- sd(values, na.rm = TRUE) / sqrt(n_samples)
   
   # Formater le résultat : moyenne (erreur standard)
   if (is.na(mean_val) || is.na(se_val)) {
     return(NA)
   } else {
     return(paste0(round(mean_val, 3), " (", round(se_val, 3), ")"))
   }
 }
 
 # Créer un nouveau dataframe pour les résultats
 result_df <- df[, c("compound", "smiles", "class", "inchikey")]
 
 # Pour chaque spagg à conserver, calculer la moyenne et SE
 for (spagg in spagg_to_keep) {
   spagg_results <- character(nrow(df))
   
   for (i in 1:nrow(df)) {
     spagg_results[i] <- calculate_mean_se(df[i, ], spagg, df$compound[i])
   }
   
   result_df[[paste0("mean_", spagg)]] <- spagg_results
 }
 
 result <- result_df |>  dplyr::filter(!dplyr::if_all(5:50, ~ is.na(.)))
 return(result)
}
