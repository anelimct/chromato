
read_excel_articles <- function(file) {
  readxl::read_xlsx(file, sheet = 2) 
}

select_iso_mono <- function(data){
  variantes_monoterpenes <- c("monoterpene", "monoterpens", "monterpene", "monterpenes", "monoterp", "monoterpeness")
  
  # Standardize Compound names and filter relevant compounds
  data <- data |> 
    dplyr::mutate(Compound = ifelse(Compound %in% variantes_monoterpenes, "monoterpenes", Compound)) |> 
    dplyr::filter(Compound == "isoprene" | Compound == "monoterpenes")
}