#' Title
#'
#' @param paradise_reports_sbtr_blanks_list, une liste de tableau content les aires des pics pour chaques molécules (lignes) par samples(colonnes), les aires sont "vraiment" lié a l'espèces car le sblancs ont été soustraits  
#' @param calib_quanti, table with ng/µg of standards in a calib, une ligne par couples calib/standard = une calib est présente plusieurs fois puisqu'elle contient plusieurs standards sauf pour isoprene 
#' @param table_calib_btw_session 
#'
#' @returns
#' @export
#'
#' @examples
area_to_quanti_iso<- function(paradise_reports_sbtr_blanks_list, calib_quanti, table_calib_btw_session) {
  # Get the list of samples to process
  samples_list <- samples_paradise(paradise_reports_sbtr_blanks_list, calib_quanti)
  
  # Initialize a list to store the resulting data frames
  quanti_list <- list()
  
  # Process each dataframe in paradise_reports_list
  for (name in names(paradise_reports_sbtr_blanks_list)) {
    paradise_report <- paradise_reports_sbtr_blanks_list[[name]]
    
    # Get the calibration slope for the current session
    calib_slope <- table_calib_btw_session[table_calib_btw_session$session == name, ]
    
    # Get the list of samples for the current dataframe
    samples <- samples_list[[name]]
    
    # Add ".CDF" to the end of each sample name
    samples <- paste0(samples, ".CDF")
    
    # Iterate over each row of the paradise_report dataframe
    for (i in 1:nrow(paradise_report)) {
      # Get the New_CAS_Number for the current row
      calib_based_on <- "78-79-5"

      
      # Find the corresponding row in calib_slope
      matching_calib <- calib_slope[calib_slope$New_CAS_Number == calib_based_on, "slope"]
      
      # Check if a matching calibration slope was found
      if (nrow(matching_calib) > 0) {
        matching_calib_value <- matching_calib[[1]]
        
        # Compute quanti from calibration for the corresponding samples
        for (sample in samples) {
          
          sample_value <- paradise_report[[sample]][i]
          
          paradise_report[[sample]][i] <- sample_value / matching_calib_value
        }
      }
    }
    
    # Store the resulting data frame in the list
    quanti_list[[name]] <- paradise_report
  }
  
  return(quanti_list)
}






area_to_quanti_mono<- function(paradise_reports_sbtr_blanks_list, calib_quanti, table_calib_btw_session) {
  # Get the list of samples to process
  samples_list <- samples_paradise(paradise_reports_sbtr_blanks_list, calib_quanti)
  
  # Initialize a list to store the resulting data frames
  quanti_list <- list()
  
  table_calib_btw_session <- dplyr::distinct(table_calib_btw_session)
  
  
  for (name in names(paradise_reports_sbtr_blanks_list)) {
    paradise_report <- paradise_reports_sbtr_blanks_list[[name]]|>  dplyr::rename(compound = New_CAS_Name) 
    
    # Get the calibration slope for the current session, one line for each standrad present
    calib_slope <- table_calib_btw_session[table_calib_btw_session$session == name, ] |>  dplyr::distinct() |> 
      dplyr::mutate(
        fingerdist = stringr::str_c("fingerDisMat.", New_CAS_Name) |>
          stringr::str_replace_all("[-|,]", ".")
      )
    
    standards <- calib_slope$fingerdist
    
    dist <- as.data.frame(chemodiv::compDis(paradise_report, type = "PubChemFingerprint" )) |> 
      dplyr::select(dplyr::all_of(standards))
    dist$fingerdist <- apply(dist, 1, function(x) {
      names(x)[which.min(x)]
    })
    
    calib_based_on <- dplyr::left_join(dist, calib_slope, by = "fingerdist" ) |>  dplyr::select(New_CAS_Number) |>  dplyr::rename(calib_based_on = New_CAS_Number )
    
    paradise_report <- cbind(paradise_report, calib_based_on) 
    # Get the list of samples for the current dataframe
    samples <- samples_list[[name]]
    
    # Add ".CDF" to the end of each sample name
    samples <- paste0(samples, ".CDF")
    
    # Iterate over each row of the paradise_report dataframe
    for (i in 1:nrow(paradise_report)) {
      # Get the New_CAS_Number for the current row
      calib_based_on <- paradise_report$calib_based_on[i]
      
      # Find the corresponding row in calib_slope
      matching_calib <- calib_slope[calib_slope$New_CAS_Number == calib_based_on, "slope"]
      
      # Check if a matching calibration slope was found
      if (nrow(matching_calib) > 0) {
        matching_calib_value <- matching_calib[[1]]
        
        # Compute quanti from calibration for the corresponding samples
        for (sample in samples) {
          
          sample_value <- paradise_report[[sample]][i]
          
          paradise_report[[sample]][i] <- sample_value / matching_calib_value
        }
      }
    }
    
    # Store the resulting data frame in the list
    quanti_list[[name]] <- paradise_report
  }
  
  return(quanti_list)
  
}




compute_ER <- function(paradise_reports_quanti_list, bvocs_samples, calib_quanti) {
  # Initialize a list to store the transformed data frames
  ER_reports_list <- list()
  
  samples_list <- samples_paradise(paradise_reports_quanti_list, calib_quanti)
  
  
  
  # Iterate over each data frame in the result list
  for (name in names(paradise_reports_quanti_list)) {
    paradise_report <- paradise_reports_quanti_list[[name]]
    
    # Get the sample columns, excluding those starting with "B_"
    samples <- samples_list[[name]]
    
    samples <- paste0(samples, ".CDF")
    samples <- samples[!grepl("^B_", samples)]
    
    
    # Add ".CDF" to the ID in bvocs_samples
    bvocs_samples$ID_CDF <- paste0(bvocs_samples$ID, ".CDF")
    
    # Iterate over each row of the result data frame
    for (i in 1:nrow(paradise_report)) {
      for (sample in samples) {
        # Find the matching row in bvocs_samples
        matching_row <- bvocs_samples[bvocs_samples$ID_CDF == sample, ]
        
        if (nrow(matching_row) > 0) {
          # Get the Leaves_DM value
          leaves_dm <- matching_row$Leaves_DM * 10^-3
          sampling_duration <- matching_row$Durée/60  # in hour
          
          
          # Transform the value using the given formula
          # quantité *60 = flow entrant en heure / biomasse * Temps échantillonnage en heure ( ) * flow sortant en heure
          paradise_report[i, sample] <- paradise_report[i, sample] * 60 / leaves_dm * sampling_duration * 6
        }
      }
    }
    
    # Store the transformed data frame in the list
    ER_reports_list[[name]] <- paradise_report
  }
  
  return(ER_reports_list)
}


mark_values <- function(er_list, paradise_reports_list, lod_3x = 193500, lod_10x = 645000) {
  marked_er <- list()
  
  for (name in names(er_list)) {
    # Obtenir les données ER et les données brutes correspondantes
    er_data <- er_list[[name]]
    raw_data <- paradise_reports_list[[name]]
    cols_to_keep <-er_data |>  dplyr::select( dplyr::any_of( c("Compound Name", "Match Quality", "Compound ID", "Est. Retention Time (min)" , "Hit 1: Probability" , "Hit 1: CAS", "compound", "calib_based_on",              
"smiles", "inchikey", "calib_based_on", "New_CAS_Name")))
    
    # Identifier les colonnes d'échantillons (.CDF)
    sample_cols <- grep("\\.CDF$", names(raw_data), value = TRUE)
    
    
  
    er_data <- er_data |>  dplyr::select(dplyr::any_of( sample_cols))
    raw_data <- raw_data |>  dplyr::select(dplyr::any_of( sample_cols))
    
     
    
    # Vérifier la correspondance des dimensions
    if (!identical(dim(er_data), dim(raw_data))) {
      stop("Les dimensions des données ER et brutes ne correspondent pas pour ", name)
    }
    
    # Créer une copie pour modification
    marked_data <- er_data
    
    # Appliquer les seuils
    for (col in sample_cols) {
      # Remplacer les valeurs ER selon les seuils des données brutes
      marked_data[[col]] <- ifelse(raw_data[[col]] < lod_3x, "nd",
                                   ifelse(raw_data[[col]] < lod_10x, "tr", 
                                          marked_data[[col]]))
      
    }
    marked_data <- cbind(marked_data, cols_to_keep)
    marked_er[[name]] <- marked_data
  }
  
  return(marked_er)
}



sum_terpenoids <- function(df) {
  # Ensure we have required columns
  if (!all(c("superclass", "superclass2") %in% colnames(df))) {
    stop("Dataframe must contain 'superclass' and 'superclass2' columns")
  }
  
  # Process the dataframe
  result <- df  |>
    # Select relevant columns
    dplyr::select(`compound name`, superclass, superclass2, matches("\\.CDF$"))  |>
    dplyr::select(-matches("^B_"), -contains("calib")) |> 
    # Pivot to long format (samples as rows)
    tidyr::pivot_longer(
      cols = -c(`compound name`, superclass, superclass2),
      names_to = "Sample",
      values_to = "Value"
    )  |>
    # Convert values to numeric (handle 'nd'/'tr')
    dplyr::mutate(
      Sample = stringr::str_to_upper(Sample),
      Value = dplyr::case_when(
        Value %in% c("nd", "tr") ~ 0,
        TRUE ~ as.numeric(Value)
      )
    )  |>
    # Classify `compound name`s
    dplyr::mutate(
      Terpenoid_Type = dplyr::case_when(
        superclass == "Monoterpenoids" | superclass2 == "Monoterpenoids" ~ "Monoterpenoids",
        superclass == "Sesquiterpenoids" | superclass2 == "Sesquiterpenoids" ~ "Sesquiterpenoids",
        TRUE ~ "Other"
      )
    )  |>
    # Filter and sum by type
    dplyr::filter(Terpenoid_Type != "Other")  |>
    dplyr::group_by(Sample, Terpenoid_Type)  |>
    dplyr::summarise(Total = sum(Value, na.rm = TRUE), .groups = "drop")  |>
    dplyr::ungroup() |> 
    dplyr::group_by(Sample) |> 
    # Pivot back to wide format
    tidyr::pivot_wider(
      names_from = Terpenoid_Type,
      values_from = Total
    )  |>
    dplyr::ungroup() 
  
  if (!"Monoterpenoids" %in% names(result)) result$Monoterpenoids <- 0
  if (!"Sesquiterpenoids" %in% names(result)) result$Sesquiterpenoids <- 0
  
  result<- result |>
    dplyr::select(Sample, Monoterpenoids, Sesquiterpenoids)
  
  return(result)
}


sum_terpenoids_across_reports <- function(reports_list) {
  # Process each dataframe in the list
  lapply(reports_list, function(df) {
    # Create superclass2 if missing
    if (!"superclass2" %in% colnames(df)) df$superclass2 <- NA
    
    # Apply your existing sum_terpenoids function
    sum_terpenoids(df)
  }) |>  
    dplyr::bind_rows(.id = "Session")  # Combine results with session names
}

sum_isoprene_across_reports <- function(reports_list) {
  # Internal conversion function
  convert_nd_tr <- function(x) {
    if (is.character(x)) {
      x[x %in% c("nd", "tr")] <- "0"
    }
    as.numeric(x)
  }
  
  # Main processing function
  process_df <- function(df) {
    # Get sample columns (ending with .CDF, excluding blanks/calibrations)
    sample_cols <- grep("\\.CDF$", names(df), value = TRUE)
    sample_cols <- sample_cols[!grepl("^B_|calib", sample_cols)]
    
    # Convert values and filter isoprene
    df %>%
      # Convert nd/tr to 0 and make numeric
      dplyr::mutate(across(all_of(sample_cols), convert_nd_tr)) %>%
      # Filter isoprene compounds
      dplyr::filter(New_CAS_Name %in% c("Isoprene", "78-79-5")) %>%
      # Select and reshape data
      dplyr::select(all_of(sample_cols)) %>%
      tidyr::pivot_longer(
        cols = everything(),
        names_to = "Sample",
        values_to = "Isoprene"
      ) %>%
      dplyr::mutate(Sample = toupper(Sample))
  }
  
  # Apply to all dataframes in list
  purrr::map_dfr(names(reports_list), ~ {
    process_df(reports_list[[.x]]) %>% 
      dplyr::mutate(Session = .x) %>%
      dplyr::select(Session, Sample, Isoprene)
  })
}


keep_terpenoids_across_reports <- function(reports_list) {
  # Process each dataframe in the list
  lapply(reports_list, function(df) {
    # Create superclass2 if missing
    if (!"superclass2" %in% colnames(df)) df$superclass2 <- NA
    
    df |> dplyr::filter( 
                  superclass == "Monoterpenoids" | superclass2 == "Monoterpenoids" |
                    superclass == "Sesquiterpenoids" | superclass2 == "Sesquiterpenoids") |> 
      dplyr::mutate(
        class = dplyr::case_when(
          (superclass == "Monoterpenoids" | superclass2 == "Monoterpenoids") & stringr::str_detect(smiles, "O") ~ "Oxygenated-monoterpenes",
          (superclass == "Sesquiterpenoids" | superclass2 == "Sesquiterpenoids") & stringr::str_detect(smiles, "O") ~ "Oxygenated-sesquiterpenes",
          
          superclass == "Monoterpenoids" | superclass2 == "Monoterpenoids" ~ "Monoterpenes",
          superclass == "Sesquiterpenoids" | superclass2 == "Sesquiterpenoids" ~ "Sesquiterpenes",
          #Garder la valeur existante si aucune condition ne correspond
          TRUE ~ class
        )
      )
      
  })
}

keep_above_5_ng <- function(reports_list, threshold = 5){
  for (name in names(reports_list)) {
    sample_cols <- grep("^[A-Za-z]{2,4}_", names(reports_list[[name]]), value = TRUE)
    
    # Create a copy for modification
    df <- reports_list[[name]]
    
    for (col in sample_cols) {
      if (is.numeric(df[[col]])) {
        # Replace values below threshold with NA
        df[[col]][df[[col]] < threshold] <- "nd"
      }
    }
    
    reports_list[[name]] <- df
  }
  return(reports_list)
}


save_ER_xlsx <- function(reports_list, suffixe_facultatif = NULL){
  
  for (name in names(reports_list)) {
    
    # Vérification sécurisée des colonnes
    expected_cols <- c("compound", "calib_based_on", "class", "est. retention time (min)")
    available_cols <- intersect(expected_cols, names(reports_list[[name]]))
    
    # Colonnes sample (doivent commencer par 4 lettres + _)
    sample_cols <- grep("^[A-Za-z]{4}_", names(reports_list[[name]]), value = TRUE)
    
    if (length(sample_cols) == 0) {
      warning("No sample columns found in ", name, ". Skipping.")
      next
    }
    
    # Sélection sécurisée des colonnes
    df <- reports_list[[name]] |> 
      dplyr::select(
        dplyr::any_of(available_cols),  # Seulement les colonnes disponibles
        dplyr::all_of(sample_cols)      # Toutes les colonnes sample doivent exister
      ) 
    
    # Vérifier si 'class' existe pour l'arrangement
    if ("class" %in% names(df)) {
      df <- df |> dplyr::arrange(class)
    }
    
    # Convertir en caractères - version sécurisée
    df <- df |> 
      dplyr::mutate(
        dplyr::across(
          .cols = dplyr::everything(),
          .fns = as.character
        )
      )
    
    # Identifier les colonnes numériques (les colonnes sample)
    numeric_cols <- sample_cols
    
    # Calcul des sommes par classe - seulement si 'class' existe
    if ("class" %in% names(df)) {
      sum_by_class <- df |>
        dplyr::mutate(
          dplyr::across(
            .cols = dplyr::all_of(numeric_cols),
            .fns = ~ dplyr::case_when(
              . %in% c("nd", "ND", "tr", "TR") ~ "0",
              TRUE ~ .
            )
          )
        ) |>
        dplyr::mutate(
          dplyr::across(
            .cols = dplyr::all_of(numeric_cols),
            .fns = ~ suppressWarnings(as.numeric(.))
          )
        ) |>
        dplyr::group_by(class) |>
        dplyr::summarise(
          dplyr::across(
            .cols = dplyr::all_of(numeric_cols),
            .fns = ~sum(., na.rm = TRUE),
            .names = "{.col}"
          ),
          compound = "SUM_ng",
          calib_based_on = NA_character_,
          `est. retention time (min)` = NA_real_,
          .groups = "drop"
        )
      
      # Réorganiser les colonnes
      if ("calib_based_on" %in% available_cols) {
        sum_by_class <- sum_by_class |> 
          dplyr::relocate(compound, calib_based_on, class, `est. retention time (min)`)
      } else {
        sum_by_class <- sum_by_class |> 
          dplyr::relocate(compound, class, `est. retention time (min)`)
      }
      
      # Conversion en µg
      sum_by_class_ug <- sum_by_class |>
        dplyr::mutate(
          dplyr::across(
            .cols = dplyr::all_of(numeric_cols),
            .fns = ~ . / 1000
          ),
          compound = "SUM_ug"
        ) 
      
      sums_combined <- dplyr::bind_rows(sum_by_class, sum_by_class_ug)
      
      # Convertir en caractères pour la fusion
      sums_combined <- sums_combined |>
        dplyr::mutate(
          dplyr::across(
            .cols = dplyr::everything(),
            .fns = as.character
          )
        )
      
      # Combiner avec les données originales
      final_df <- dplyr::bind_rows(df, sums_combined)
    } else {
      # Si pas de colonne 'class', utiliser les données originales
      final_df <- df
      warning("No 'class' column in ", name, ". Skipping summary calculations.")
    }
    
    # Nom du fichier
    if (!is.null(suffixe_facultatif)) {
      file_name <- paste0("ER_", stringr::str_extract(name, "^[^.]*"), "_", suffixe_facultatif, ".xlsx")
    } else {
      file_name <- paste0("ER_", name)
    }
    
    file_path <- file.path(here::here("outputs", "ER_reports"), file_name)
    openxlsx::write.xlsx(final_df, file = file_path)
  }
  
  return(reports_list)
}




transformer_df <- function(df) {
  # Mettre tous les noms de colonnes en minuscules
  names(df) <- tolower(names(df))
  
  # Renommer 'new_cas_name' en 'compound'
  if("new_cas_name" %in% names(df)) {
    names(df)[names(df) == "new_cas_name"] <- "compound"
  }
  
  # Ajouter la colonne 'class' avec la valeur de 'compound'
  if("compound" %in% names(df)) {
    df$class <- df$compound
  }
  
  return(df)
}