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
  list_standardized_tables <- list(compounds_table_standardized,compounds_table_standardized_mono_LetT, compounds_table_standardized_mono_T, compounds_table_standardized_mono_T_bourtsou )
  return(list_standardized_tables )
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
    ) |>  dplyr::filter(detection_rate >= 0.5) |>  dplyr::mutate(spagg = stringr::str_to_lower(spagg))
}



is_2025_sample <- function(col_name) {
  
  grepl("2025",col_name)
}


compound_mean_sp <- function(compounds_table, valid_samples_mono, times_compound_sp, include_se = TRUE){
  
  cols_to_keep <- valid_samples_mono$ID |> stringr::str_to_lower() |> stringr::str_c(".cdf", sep = "") |>
    (\(x) x[grepl("^[a-zA-Z]{4}_", x)])()
  
  df <- compounds_table |> dplyr::select(any_of(c("compound", "smiles", "class", "inchikey", paste0(cols_to_keep))))
  
  spagg_to_keep <- unique(times_compound_sp$spagg)
  
  # Helper function to check if sample is from 2025
  is_2025_sample <- function(sample_name) {
    # Adjust this logic based on how 2025 samples are identified in your data
    grepl("2025", sample_name) || grepl("25_", sample_name)
  }
  
  # Fonction pour calculer la moyenne et l'erreur standard pour un spagg donné
  calculate_mean_se <- function(compound_row, spagg, compound_name, include_se, times_compound_sp) {
    # Trouver les colonnes correspondant au spagg
    spagg_cols <- names(compound_row)[grepl(paste0("^", spagg), names(compound_row))]
    
    # If no columns found, return NA
    if (length(spagg_cols) == 0) {
      return(NA)
    }
    
    values <- as.numeric(compound_row[spagg_cols])
    
    compound_class <- compound_row$class
    
    # Trouver le nombre d'échantillons pour la division
    if (compound_class %in% c("Sesquiterpenes", "Oxygenated-sesquiterpenes")) {
      # Pour ces classes, utiliser seulement les échantillons de 2025
      sample_cols_2025 <- spagg_cols[sapply(spagg_cols, is_2025_sample)]
      n_samples <- length(sample_cols_2025)  
    } else {
      # Pour les autres classes, utiliser total_samples_species de times_compound_sp
      species_compound_info <- times_compound_sp |> 
        dplyr::filter(spagg == !!spagg, compound == !!compound_name)
      
      # Check if we found matching info
      if (nrow(species_compound_info) == 0) {
        n_samples <- length(spagg_cols)  # fallback to number of columns
      } else {
        n_samples <- species_compound_info$total_samples_species[1]
      }
    }
    
    # Ensure n_samples is valid
    if (is.na(n_samples) || n_samples == 0) {
      return(NA)
    }
    
    # Remove NA values before calculation
    valid_values <- values[!is.na(values)]
    
    if (length(valid_values) == 0) {
      return(NA)
    }
    
    # Calculer la moyenne
    mean_val <- sum(valid_values, na.rm = TRUE) / n_samples
    
    # Calculer l'erreur standard
    se_val <- sd(valid_values, na.rm = TRUE) / sqrt(n_samples)
    
    # Formater le résultat selon l'option include_se
    if (is.na(mean_val) || is.na(se_val)) {
      return(NA)
    } else {
      if (include_se) {
        return(paste0(round(mean_val, 3), " (", round(se_val, 3), ")"))
      } else {
        return(round(mean_val, 3))
      }
    }
  }
  
  # Créer un nouveau dataframe pour les résultats
  result_df <- df[, c("compound", "smiles", "class", "inchikey")]
  
  # Pour chaque spagg à conserver, calculer la moyenne et SE
  for (spagg in spagg_to_keep) {
    spagg_results <- character(nrow(df))
    
    for (i in 1:nrow(df)) {
      spagg_results[i] <- calculate_mean_se(df[i, ], spagg, df$compound[i], include_se, times_compound_sp)
    }
    
    result_df[[paste0("mean_", spagg)]] <- spagg_results
  }
  
  # FIXED: Use dynamic column selection instead of hard-coded indices
  mean_cols <- grep("^mean_", names(result_df), value = TRUE)
  
  if (length(mean_cols) > 0) {
    result <- result_df |> 
      dplyr::filter(!dplyr::if_all(dplyr::all_of(mean_cols), ~ is.na(.)))
  } else {
    result <- result_df
  }
  
  return(result)
}

compounds_samples_spagg_to_keep <- function(compounds_table, valid_samples_mono, times_compound_sp) {
  
  # Récupérer les spagg à garder
  spagg_to_keep <- unique(times_compound_sp$spagg)
  
  # Créer un pattern regex pour les spagg à garder
  spagg_pattern <- paste0("^(", paste(spagg_to_keep, collapse = "|"), ")_")
  
  # Filtrer les échantillons valides qui commencent par un spagg_to_keep
  filtered_sample_ids <- valid_samples_mono$ID |> 
    stringr::str_to_lower() |> 
    stringr::str_c(".cdf") |>
    stringr::str_subset(spagg_pattern)
  
  # Colonnes à garder
  cols_to_keep <- c("compound", "smiles", "class", "inchikey", filtered_sample_ids)
  
  # Filtrer le dataframe
  result_df <- compounds_table |> 
    dplyr::select(dplyr::any_of(cols_to_keep))
  
  # Enlever les lignes avec seulement des NAs dans les colonnes d'échantillons
  sample_cols <- setdiff(names(result_df), c("compound", "smiles", "class", "inchikey"))
  
  if (length(sample_cols) > 0) {
    result_df <- result_df |> 
      dplyr::filter(dplyr::if_any(dplyr::all_of(sample_cols), ~ !is.na(.)))
  }
  
  return(result_df)
}

## put 0 for a compound in a sample if this compound is a singelton = not found in at least 50% of samples

compounds_tabled_zeroed_singleton <- function(compounds_table, times_compound_sp) {
  # Extraire les couples espèce-composé valides
  valid_compounds_sp <- times_compound_sp |> 
    dplyr::select(spagg, compound) |> 
    dplyr::distinct()
  
  # Identifier les colonnes d'échantillons
  sample_cols <- names(compounds_table)[grepl("^[a-zA-Z]{4}_", names(compounds_table))]
  
  # Appliquer le filtrage pour chaque colonne échantillon
  filtered_table <- compounds_table
  
  for(sample_col in sample_cols) {
    species_code <- substr(sample_col, 1, 4)
    
    valid_compounds <- valid_compounds_sp |> 
      dplyr::filter(spagg == species_code) |> 
      dplyr::pull(compound)
    
    if(length(valid_compounds) > 0) {
      filtered_table <- filtered_table |> 
        dplyr::mutate(
          !!sample_col := ifelse(compound %in% valid_compounds, 
                                 .data[[sample_col]],  # Garder la valeur d'origine
                                 NA)                   # Mettre NA si non valide
        )
    } else {
      filtered_table[[sample_col]] <- NA
    }
  }
  
  return(filtered_table)
}




# Fonction pour calculer les sommes par classe avec écart-type
calculate_class_sums <- function(result_df, include_se = TRUE) {
  
  # Fonction pour extraire la moyenne et le SE des chaînes formatées
  extract_mean_se <- function(x) {
    if (is.na(x) || x == "NA") {
      return(list(mean = NA, se = NA))
    }
    # Extraire la moyenne (premier nombre avant l'espace et parenthèse)
    mean_val <- as.numeric(stringr::str_extract(x, "^[-+]?[0-9]*\\.?[0-9]+"))
    
    # Extraire l'erreur standard (nombre entre parenthèses)
    se_match <- stringr::str_match(x, "\\(([-+]?[0-9]*\\.?[0-9]+)\\)")
    se_val <- ifelse(is.na(se_match[1,2]), NA, as.numeric(se_match[1,2]))
    
    return(list(mean = ifelse(is.na(mean_val), 0, mean_val), 
                se = se_val))
  }
  
  # Identifier les colonnes de moyenne
  mean_cols <- grep("^mean_", names(result_df), value = TRUE)
  
  # Créer un dataframe temporaire pour les calculs
  temp_df <- result_df[, c("class", mean_cols)]
  
  # Extraire les valeurs numériques pour chaque colonne
  for (col in mean_cols) {
    extracted_vals <- lapply(temp_df[[col]], extract_mean_se)
    temp_df[[paste0("mean_", col)]] <- sapply(extracted_vals, function(x) x$mean)
    temp_df[[paste0("se_", col)]] <- sapply(extracted_vals, function(x) x$se)
  }
  
  # Colonnes extraites pour les moyennes et SE
  mean_extracted_cols <- grep("^mean_mean_", names(temp_df), value = TRUE)
  se_extracted_cols <- grep("^se_mean_", names(temp_df), value = TRUE)
  
  # Calculer la somme des moyennes par classe
  class_sums <- temp_df %>%
    dplyr::group_by(class) %>%
    dplyr::summarise(
      dplyr::across(
        all_of(mean_extracted_cols),
        ~sum(., na.rm = TRUE),
        .names = "sum_{.col}"
      )
    )
  
  # Calculer l'erreur standard de la somme (propagation d'erreur)
  class_se_sums <- temp_df %>%
    dplyr::group_by(class) %>%
    dplyr::summarise(
      dplyr::across(
        all_of(se_extracted_cols),
        ~sqrt(sum(.^2, na.rm = TRUE)),  # Propagation d'erreur pour la somme
        .names = "se_{.col}"
      )
    )
  
  # Fusionner les résultats
  class_results <- class_sums %>%
    dplyr::left_join(class_se_sums, by = "class")
  
  # Formater les résultats
  result_class <- data.frame(class = class_results$class)
  
  # Extraire les noms des espèces à partir des noms de colonnes
  species_names <- stringr::str_remove(mean_cols, "^mean_")
  
  for (species in species_names) {
    sum_col <- paste0("sum_mean_mean_", species)
    se_col <- paste0("se_se_mean_", species)
    
    sums <- class_results[[sum_col]]
    ses <- class_results[[se_col]]
    
    formatted_results <- character(length(sums))
    
    for (i in seq_along(sums)) {
      if (is.na(sums[i]) || sums[i] == 0) {
        formatted_results[i] <- NA
      } else {
        if (include_se) {
          formatted_results[i] <- paste0(round(sums[i], 3), " (", round(ses[i], 3), ")")
        } else {
          formatted_results[i] <- round(sums[i], 3)
        }
      }
    }
    
    result_class[[paste0("sum_", species)]] <- formatted_results
  }
  
  
  
  dir_path <- here::here("outputs", "ER_means")
  if (!dir.exists(dir_path)) {
    dir.create(dir_path, recursive = TRUE)
  }
  
  # Sauvegarder le tableau
  file_path <- file.path(dir_path, "means_by_class.xlsx")
  openxlsx::write.xlsx(result_class, file = file_path)
  
  return(result_class)
  
  
}
