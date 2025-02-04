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

update_lib <- function (new_lib, old_correspondant_table_with_RI) {
  old <- unlist(old_correspondant_table_with_RI$Merged_CAS_Numbers)
  new_lib <- new_lib |> dplyr::filter(!`CAS Registry Number` %in% old) #enlever toutes les lignes qui sont déja dans old
  
  updated_table <- old_correspondant_table_with_RI %>%
    left_join(new_lib, by = "Canonical SMILES") %>%
    rowwise() %>%
    mutate(
      Merged_CAS_Numbers = ifelse(
        !is.na(`CAS Registry Number.y`),
        paste(unique(c(unlist(strsplit(Merged_CAS_Numbers, ", ")), `CAS Registry Number.y`)), collapse = ", "),
        Merged_CAS_Numbers
      )
    ) %>%
    select(-`CAS Registry Number.y`) %>%
    rename(`CAS Registry Number` = `CAS Registry Number.x`)
  
  
  #case when same canonical smiles add cas num to Merged_CAS_numbers, and names to merged names
}
