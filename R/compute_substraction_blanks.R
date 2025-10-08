subtract_blanks_from_samples <- function(paradise_reports_list, bvocs_samples, calib_files) {
  
  # Obtenir la liste des samples propres avec votre fonction
  samples_list <- samples_paradise(paradise_reports_list, calib_files)
  
  # Fonction pour trouver les blanks correspondants à un sample
  find_matching_blanks <- function(sample_id, bvocs_samples) {
    # Extraire les infos du sample
    sample_row <- bvocs_samples |> 
      dplyr::filter(ID == sample_id)
    
    if (nrow(sample_row) == 0) {
      return(character(0))
    }
    
    # Trouver les blanks correspondants (même date, caisse, zone géographique)
    matching_blanks <- bvocs_samples |> 
      dplyr::filter(
        Taxon == "BLANC",
        Date == sample_row$Date,
        N..caisse == sample_row$N..caisse,
        abs(Latitude..WGS84. - sample_row$Latitude..WGS84.) < 0.03,
        abs(Longitude..WGS84. - sample_row$Longitude..WGS84.) < 0.03
      ) |> 
      dplyr::pull(ID)
    
    return(matching_blanks)
  }
  
  # Traiter chaque batch
  results <- list()
  
  for (batch_name in names(paradise_reports_list)) {
    paradise_report <- paradise_reports_list[[batch_name]]
    batch_samples <- samples_list[[batch_name]]
    
    # FILTRE IMPORTANT : exclure les blanks de la liste des samples à traiter
    real_samples <- batch_samples[!grepl("^B_", batch_samples)]
    
    # Créer une copie pour les résultats
    result_df <- paradise_report
    
    # Pour chaque VRAI sample du batch (pas les blanks)
    for (sample_id in real_samples) {
      sample_col <- paste0(sample_id, ".CDF")
      
      # Vérifier que la colonne existe dans le rapport
      if (!sample_col %in% colnames(paradise_report)) {
        warning(paste("Colonne", sample_col, "non trouvée dans", batch_name))
        next
      }
      
      # Étape 1: Trouver les blanks correspondants dans les métadonnées
      matching_blanks <- find_matching_blanks(sample_id, bvocs_samples)
      
      if (length(matching_blanks) == 0) {
        warning(paste("Aucun blank trouvé pour le sample", sample_id, "dans", batch_name))
        next
      }
      
      # Étape 2: Vérifier que les blanks existent dans le rapport
      blank_cols <- paste0(matching_blanks, ".CDF")
      available_blanks <- blank_cols[blank_cols %in% colnames(paradise_report)]
      
      if (length(available_blanks) == 0) {
        warning(paste("Blanks trouvés mais absents du rapport pour", sample_id, "dans", batch_name))
        next
      }
      
      # Soustraire la moyenne des blanks pour chaque composé
      for (i in 1:nrow(paradise_report)) {
        blank_values <- as.numeric(paradise_report[i, available_blanks])
        blank_mean <- mean(blank_values, na.rm = TRUE)
        
        if (!is.na(blank_mean)) {
          original_value <- as.numeric(paradise_report[i, sample_col])
          corrected_value <- original_value - blank_mean
          result_df[i, sample_col] <- max(corrected_value, 0)  # Éviter les valeurs négatives
        }
      }
      
      # Message de confirmation
      message(paste("Soustraction appliquée pour", sample_id, "avec", length(available_blanks), "blanks"))
    }
    
    results[[batch_name]] <- result_df
  }
  
  return(results)
}
