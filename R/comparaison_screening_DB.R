

comparison_screening_DB <- function(sp_screening, db) {
  
  
conditions <- db |> dplyr::filter(
  gragg %in% sp_screening$gragg) |>  dplyr::filter(Compound %in% c("monoterpenes", "isoprene", "sesquiterpenes", "sesquiterepenes") ) |> 
  dplyr::filter(Emission_unit_comp == "microg" & Emission_unit_leaf == "g" &
                  Emission_unit_time == "h") |>  dplyr::filter(Standardized %in% c("true", "TRUE") | Temperature == "30" & Temperature_unit == "C"| Temperature == "303" & Temperature_unit == "K"
                                                 | Temperature == "303.15" & Temperature_unit == "K")
}
