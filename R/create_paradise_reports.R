

read_paradise_4 <- function(files_path, CAS_library, option = "") {
  # Initialize the list to store data frames
  list_pardise_files <- list()
  
  # Iterate over each file path to have the full path
  for (i in seq_along(files_path)) {
    file_path <- paste0(here::here("data", "paradise_reports"), "/", option,  files_path[i])
    
    # Read the Excel sheets
    peak_area <- readxl::read_xlsx(file_path, sheet = 2)
    top_nist <- readxl::read_xlsx(file_path, sheet = 4) |>
      dplyr::select("Compound Name","Match Quality",               "Compound ID",
                    "Interval ID",                 "Model Order",                 "Component #",
                    "Est. Retention Time (min)"   ,  "Hit 1: Probability", "Hit 1: CAS")
    
    # Join the data frames
    data <- dplyr::left_join(peak_area, top_nist, by = c("Compound Name", "Match Quality", "Compound ID", "Interval ID", "Model Order","Component #", "Est. Retention Time (min)"))
    

    # Filter CAS_library and extract relevant information
    result <- data |>
      dplyr::rowwise() |> 
      dplyr::mutate(
        match = list(CAS_library |>
                       dplyr::filter(sapply(Merged_CAS_Numbers, function(x) `Hit 1: CAS` %in% x)))
      ) |> #creates a new column named match in the data frame. For each row, it filters the CAS_library data frame based on whether the Hit 1: CAS value is found in the Merged_CAS_Numbers list
      dplyr::ungroup() |>
      dplyr::mutate(
        New_CAS_Name = sapply(match, function(x) if (nrow(x) > 0) x$New_CAS_Name[1] else NA),
        RI = sapply(match, function(x) if (nrow(x) > 0) x$RI[1] else NA),
        nb_C_n = sapply(match, function(x) if (nrow(x) > 0) x$nb_C_n[1] else NA),
        nb_C_n_s = sapply(match, function(x) if (nrow(x) > 0) x$nb_C_n_s[1] else NA), 
        smiles = sapply(match, function(x) if (nrow(x) > 0) x$`Canonical SMILES`[1] else NA), 
        inchikey = sapply(match, function(x) if (nrow(x) > 0) x$inchikey[1] else NA) 
      ) |> #extract info in CAS lib for the matching row, prendre le premier match sinon NA
      dplyr::select(-match)

    # Store the result in the list
    list_pardise_files[[basename(files_path[i])]] <- result
    
  }
  
  return(list_pardise_files)
}



#' Title
#'
#' @param pradise_list, list containing data frames with paradise reports 
#' @param calib_files , file with the firt colum being names of calib files
#'
#' @returns Cette fonction prendre la liste des paradise reports et retourn le nom des tubes présents dans le report (sans.CDF ni .D)
#' @export
#'
#' @examples
samples_paradise <- function(paradise_reports_list, calib_files){
  # dans le paradise report il y a une colonne pour chaque samples et donc il y a aussi des colonnes pour les calib que l'on veut ignorées
  calib_files <- stringr::str_replace(calib_files[,1], ".D", ".CDF")
  columns_to_exclude <- c("Compound Name", "New_CAS_Name", "Match Quality", "Compound ID", "Interval ID", "Model Order" , "Component #", "Est. Retention Time (min)", "Comments", "Hit 1: Probability" , "Hit 1: CAS", calib_files, "RI",  "nb_C_n", "nb_C_n_s", "calib_based_on", "smiles", "inchikey" )
  
  # Initialiser une liste pour stocker les noms de samples présents dans un paradise batch
  paradise_samples_list <- list()
  
  # Parcourir chaque DataFrame dans la liste
  for (name in names(paradise_reports_list)) {
    data <- paradise_reports_list[[name]]
    paradise_samples <- setdiff(colnames(data), columns_to_exclude)
    paradise_samples <- stringr::str_remove(paradise_samples, ".CDF")
    paradise_samples_list[[name]] <- paradise_samples
  }
  return(paradise_samples_list)
}

