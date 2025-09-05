area_to_quanti_iso <- function(paradise_reports_sbtr_blanks_list, calib_quanti, table_calib_btw_session) {
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
          paradise_report[i, sample] <- paradise_report[i, sample] * 60 / leaves_dm * 0.25 * 6
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
    
    # Identifier les colonnes d'échantillons (.CDF)
    sample_cols <- grep("\\.CDF$", names(raw_data), value = TRUE)
    
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
  }) %>% 
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
