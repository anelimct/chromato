create_list_dataframes_blanks <- function(bvocs_samples){
  unique(bvocs_samples$batch)
  
  list_blanks_tables <- list()
  
  bvocs_samples<-  bvocs_samples |> dplyr::filter(!batch %in% c(NA, "wNA_NA", "NA"))  |> 
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
  
  cols_to_include <- c( "Compound Name","Match Quality",  "Est. Retention Time (min)", "New_CAS_Name", "calib_based_on", "smiles")
  
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



