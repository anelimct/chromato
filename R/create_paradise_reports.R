read_paradise <- function(files_path, CAS_library){
  
  for (i in seq_along(files_path)) {
    
    file_path <- files_path[i]
    
    
    peak_area <- readxl::read_xlsx( file_path, sheet = 2)
    top_nist <- readxl::read_xlsx( file_path, sheet = 4) |> 
      dplyr::select("Compound Name" , "Hit 1: Probability", "Hit 1: CAS" )
    
    data <- dplyr::left_join(peak_area, top_nist)
    
    result <- data %>%
      rowwise() %>%
      # Filtrer CAS_library pour chaque ligne de data
      mutate(match = list(CAS_library %>%
                            filter(sapply(`Merged_CAS_Numbers`, function(x) `Hit 1: CAS` %in% x)))) %>%
      ungroup() %>%
      # Extraire les colonnes souhaitées du match
      mutate(
        RI = ifelse(nrow(match[[1]]) > 0, match[[1]]$RI[[1]], NA),
        nb_C_n = ifelse(nrow(match[[1]]) > 0, match[[1]]$nb_C_n[[1]], NA), 
        nb_C_n_s = ifelse(nrow(match[[1]]) > 0, match[[1]]$nb_C_n_s[[1]], NA)
      ) %>%
      select(-match)  
    
    
  }
}
