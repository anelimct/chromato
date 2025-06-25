create_list_dataframes_blanks <- function(bvocs_samples){
  unique(bvocs_samples$batch)
  
  list_blanks_tables <- list()
  
  bvocs_samples<-  bvocs_samples |> dplyr::filter(!batch %in% c(NA, "NA"))  |> 
    dplyr::group_by(batch) |>dplyr::rename(list_ID_blank = `list(ID_blank)`) |>
    dplyr::mutate(group_blanks = sapply(list_ID_blank, function(x) {
      unique(sapply(unlist(x), function(y) {
        if (stringr::str_detect(y, "H")) {
          stringr::str_sub(y, 12, 19)
        } else {
          stringr::str_sub(y, 10, 17)
        }
      }))
    })) |>  dplyr::ungroup()
  
  for (batch_x in unique(bvocs_samples$batch)) {

      # Filter the data for the current batch
    batch_data <- bvocs_samples |>  dplyr::filter(batch == batch_x)

      # Get unique group_blanks for the current batch
      unique_groups <- unique(batch_data$group_blanks)

      # Initialize an empty dataframe with NA values
      s <- data.frame(matrix(NA, ncol = length(unique_groups), nrow = 4))

      # Set the column names to be the unique group_blanks
      colnames(s) <- unique_groups

      for (i in 1 : length( unique_groups)) {
        # Get the list_ID_blank values for the current group
        id_blanks <- unique(unlist(batch_data$list_ID_blank[batch_data$group_blanks == paste0(unique_groups[i])]))
        id_blanks<- paste0(id_blanks, ".CDF")
        
        
        s[1:4, unique_groups[i]] <- c(id_blanks, rep(NA, 4 - length(id_blanks)))
      }

      list_blanks_tables[[paste0(batch_x, ".xlsx")]] <- s
  }
  return(list_blanks_tables)
}

sort_list_blanks_tables <- function(list_blanks_tables) {
  # Helper function to keep only unique columns in a dataframe
  keep_unique_columns <- function(df) {
    df[, !duplicated(colnames(df))]
  }
  
  # Combine w11_2023.xlsx and w12_2023.xlsx, and name it c15_2023.xlsx
  if ("w11_2024.xlsx" %in% names(list_blanks_tables) & "w12_2024.xlsx" %in% names(list_blanks_tables)) {
    combined_df <- cbind(list_blanks_tables[["w11_2024.xlsx"]], list_blanks_tables[["w12_2024.xlsx"]])
    list_blanks_tables[["c5_2023.xlsx"]] <- keep_unique_columns(combined_df)
  }
  
  # Remove w11_2023.xlsx and w12_2023.xlsx from the list
  list_blanks_tables[["w11_2024.xlsx"]] <- NULL
  list_blanks_tables[["w12_2024.xlsx"]] <- NULL
  
  # Remove wNA_NA.xlsx from the list
  list_blanks_tables[["wNA_NA.xlsx"]] <- NULL
  
  # Combine w27_2023.xlsx, w28_2023.xlsx, and w30_2023.xlsx, and name it c5_2023.xlsx
  if ("w27_2023.xlsx" %in% names(list_blanks_tables) & "w28_2023.xlsx" %in% names(list_blanks_tables) & "w30_2023.xlsx" %in% names(list_blanks_tables)) {
    combined_df <- cbind(list_blanks_tables[["w27_2023.xlsx"]], list_blanks_tables[["w28_2023.xlsx"]], list_blanks_tables[["w30_2023.xlsx"]])
    list_blanks_tables[["c15_2023.xlsx"]] <- keep_unique_columns(combined_df)
  }
  
  # Remove w27_2023.xlsx, w28_2023.xlsx, and w30_2023.xlsx from the list
  list_blanks_tables[["w27_2023.xlsx"]] <- NULL
  list_blanks_tables[["w28_2023.xlsx"]] <- NULL
  list_blanks_tables[["w30_2023.xlsx"]] <- NULL
  
  return(list_blanks_tables)
}


paradise_reports_grouped_blanks <- function(paradise_reports_list, blanks_tables_list, iso = FALSE ){
  
  if (length(blanks_tables_list) != length(paradise_reports_list)) {
    stop("There is a mismatch between blanks_tables_list and paradise_reports_list.")
  }
  if (iso) {
    names(blanks_tables_list) <- sub("\\.xlsx$", "_iso.xlsx", names(blanks_tables_list))
  }
  
  list_blanks_grp_batch <- list()
  
  cols_to_include <- c( "Compound Name","Match Quality",  "Est. Retention Time (min)", "New_CAS_Name", "calib_based_on")
  
  for (name in names(blanks_tables_list)){
    table_blanks <- blanks_tables_list[[name]]
    paradise_report_session <- paradise_reports_list[[name]]
    
    s <- paradise_report_session |>  dplyr::select(cols_to_include)
    
    
    # Debugging: Print column names to verify they match
    print(paste("Processing:", name))
      
    for (group_blanks in  colnames(table_blanks)){
      blanks <- na.omit(table_blanks[,  group_blanks])
      
      blanks <- dplyr::setdiff(blanks, c("B_198952_15062023.CDF", "B_809115_15062023.CDF", "B_809140_15062023.CDF", "B_809118_28062023.CDF", "B_809101_29062023.CDF", "B_809102_06072023.CDF" , "B_809108_06072023.CDF", "B_809116_07072023.CDF" ,"B_809124_07072023.CDF", "B_809133_19072023.CDF", "B_809110_20072023.CDF" ,  "B_809123_20072023.CDF"))
      
      
      mean_col_name <- paste0(group_blanks, "_µ")
      diff_col_name <- paste0(group_blanks, "_Δ")
      
      # Create a data frame with the new columns
      data_g <- paradise_report_session |> 
        dplyr::rowwise() |> 
        dplyr::mutate(!!mean_col_name := mean(dplyr::c_across(all_of(blanks)), na.rm = TRUE),
               !!diff_col_name := max(dplyr::c_across(all_of(blanks)), na.rm = TRUE) - min(dplyr::c_across(all_of(blanks)), na.rm = TRUE)) |> 
        dplyr::ungroup() |>  dplyr::select(c(cols_to_include, mean_col_name, diff_col_name))
      
      print(colnames(data_g))
      print(length(rownames(data_g)))
            
      if (length(rownames(data_g)) > length(rownames(paradise_report_session))) {
        stop("There is too much ranks in data_g")

      }
      
      s <- dplyr::inner_join(s, data_g, by = cols_to_include)
      
    }
    list_blanks_grp_batch[[name]] <- s
  }
 return(list_blanks_grp_batch)
}



subtract_blanks_from_samples <- function(paradise_reports_list, grouped_blanks_list, calib_quanti) {
  # Get the list of samples to process
  samples_list <- samples_paradise(paradise_reports_list, calib_quanti)
  
  # Process each dataframe in paradise_reports_list
  for (name in names(paradise_reports_list)) {
    paradise_report <- paradise_reports_list[[name]]
    grouped_blanks <- grouped_blanks_list[[name]]
    
    # Get the list of samples for the current dataframe
    samples <- samples_list[[name]]
    
    # Filter out samples that start with "B_" and add ".CDF" to the end
    samples <- paste0(samples, ".CDF")
    samples <- samples[!grepl("^B_", samples)]
    
    # Iterate over each row of the paradise_report dataframe
    for (i in 1:nrow(paradise_report)) {
      # Get the New_CAS_Name for the current row
      cas_name <- paradise_report$New_CAS_Name[i]
      
      # Find the corresponding row in grouped_blanks
      matching_row <- which(grouped_blanks$New_CAS_Name == cas_name)
      
      if (length(matching_row) > 0) {
        matching_row <- matching_row[1]
        
        # Subtract the mean blank values from the corresponding samples
        for (sample in samples) {
          # Extract the date from the sample name
          date <- sub("^[^_]*_[^_]*_([^_]*)\\.CDF$", "\\1", sample)
          
          # Get the corresponding mean blank value
          mean_col_name <- paste0(date, "_µ")
          if (mean_col_name %in% colnames(grouped_blanks)) {
            mean_blank_value <- grouped_blanks[matching_row, mean_col_name]
            
            # Subtract the mean blank value from the sample column
            if (sample %in% colnames(paradise_report)) {
              paradise_report[i, sample] <- paradise_report[i, sample] - mean_blank_value
            }
          }
        }
      }
    }
    
    # Update the paradise_reports_list with the modified dataframe
    paradise_reports_list[[name]] <- paradise_report
  }
  
  return(paradise_reports_list)
}




area_to_quanti <- function(paradise_reports_sbtr_blanks_list, calib_quanti, table_calib_btw_session) {
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

