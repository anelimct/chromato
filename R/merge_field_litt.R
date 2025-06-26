merge_datasets <- function(df_ERfield_mono, df_ERfield_iso, bvocs_samples){
  
  df_iso_mono <- merge(df_ERfield_mono, df_ERfield_iso, by = "Sample" ) |> dplyr::mutate(ID = stringr::str_replace(Sample, ".CDF", "")) |> dplyr::select(ID, Monoterpenoids, Sesquiterpenoids, Isoprene) |> dplyr::rename(isoprene = Isoprene, monoterpenes = Monoterpenoids, sesquiterpenes = Sesquiterpenoids)
  
  
  
  cols_PAR <- c("PAR.début.1",        "PAR.début.2",       "PAR.milieu.1",     
               "PAR.milieu.2",  "PAR.milieu.3", "PAR.milieu.4",  "PAR.milieu.5", "PAR.milieu.6", "PAR.fin.1",     "PAR.fin.2")
  
  
  bvocs_samples_ER <- merge(bvocs_samples, df_iso_mono, by = "ID") |> dplyr::filter(is.na(Leaves_DM)) |>   dplyr::mutate(across(all_of(cols_PAR), as.numeric)) |>  dplyr::mutate(PAR_algo = rowMeans(dplyr::across( all_of(cols_PAR)), na.rm = TRUE)) |> 
    
    dplyr::mutate(mean_T = purrr::map_dbl(values_T_in, ~ mean(.x))) |> dplyr::mutate(T_algo_K = mean_T + 273,15) |> 
    dplyr::select(ID, monoterpenes, isoprene, Taxon, PAR_algo, T_algo_K ) |>  tidyr::pivot_longer(cols = c(isoprene, monoterpenes),names_to = "Compound",values_to = "Emission") |> dplyr::mutate(Standardized = FALSE) |>  apply_standardization()
  
  
}
  