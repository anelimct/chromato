#' Title
#'
#' @param files_path, list de chemins ou se trouvent les fichiers CAS téléchargés 
#'
#' @returns
#' @export
#'
#' @examples
create_library <- function(files_path){
  CAS_lib <- c()
  
  for (i in seq_along(files_path)) {
    data <- openxlsx::read.xlsx(files_path[i])
    
    
    data <- data[-1, ]  
    colnames(data) <- data[1, ]  #
    data <- data[-1, ]
    
    CAS_lib <- rbind(CAS_lib, data)
    CAS_lib <- unique(CAS_lib)
    
  }
  return(CAS_lib)
}


update_lib <- function (new_lib) {
  
  old_correspondant_table_with_RI <- readxl::read_xlsx(here::here("data", "library_cas_ri.xlsx"), sheet = 1)
  
  old <- unlist(old_correspondant_table_with_RI$Merged_CAS_Numbers)
  new_lib <- new_lib |> dplyr::filter(!`CAS Registry Number` %in% old) #enlever toutes les lignes qui sont déja dans old
  
  # permet de merged tous les noms de composants et leur numéros CAS pour ceuw qui avit déjà un Canonical smiles dans la table des correspondances
  updated_table <- old_correspondant_table_with_RI |> 
    dplyr::left_join(new_lib, by = "Canonical SMILES") |> 
    dplyr::group_by(`Canonical SMILES`) |> 
    dplyr::summarise(
           New_CAS_Number = unique(`New_CAS_Number`),
           New_CAS_Name = unique(`New_CAS_Name`),
           Merged_CAS_Numbers = list(unique(purrr::discard(c(unlist(`Merged_CAS_Numbers`), `CAS Registry Number`), is.na))),
           Merged_CAS_Names = list(unique(purrr::discard(c(unlist(`Merged_CAS_Names`), `CAS Index Name`), is.na))), 
           RI = unique(`RI`)
           )
  included <- unlist(updated_table$Merged_CAS_Numbers)
  new_lib <- new_lib |> dplyr::filter(!`CAS Registry Number` %in% included) #garder les lignes qui n'avaient pas de canonical smiles dans la correspondant table et qui n'ont donc pas été inclus dans le left join
  
  new_smiles <-  new_lib  %>%
    dplyr::group_by(`Canonical SMILES`) %>%
    dplyr::summarise(
      New_CAS_Number = dplyr::first(`CAS Registry Number`),  # Use the first CAS number in the group
      New_CAS_Name = dplyr::first(`CAS Index Name`),        # Use the first CAS name in the group
      Merged_CAS_Numbers = list(`CAS Registry Number`), Merged_CAS_Names = list(`CAS Index Name`), # Track merged CAS numbers
      .groups = "drop"
    ) |> tibble::add_column(RI = as.double(NA) )
  
  output <- dplyr::bind_rows(new_smiles, updated_table)
  writexl::write_xlsx( output, path = here::here("data", "library_cas_ri.xlsx"))
 
  return(output)
}
