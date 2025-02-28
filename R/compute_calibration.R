compare_calib_btw_reports <- function (calib_files, paradise_reports_list, paradise_reports_calib_list){
  #puisque les samples sont en .CDF et non D les renommer
  calib_files <- calib_files |> dplyr::mutate(new_name =stringr::str_replace(new_name, ".D", ".CDF") )
   
  calib_samples_mono <- unique(calib_files[stringr::str_detect( calib_files$new_name, "mono"), 1])

  cols_to_include_mono <- c(calib_samples_mono,  "New_CAS_Name", "Hit 1: CAS" )
  
  
  calib_btw_batch <- paradise_reports_list[[1]] |> dplyr::select("New_CAS_Name", "Hit 1: CAS") |> tibble::add_column(calib = as.character(NA)) 
  
  calib_btw_batch_mono <- calib_btw_batch[0,]
  
  
  calib_solo_mono <- paradise_reports_calib_list[["calib_mono.xlsx"]] |> dplyr::filter(New_CAS_Name %in% c("(±)-α-Pinene", "β-Pinene", "Limonene", "p-Cymene (8CI)" ))|>
    dplyr::select(tidyselect::all_of(cols_to_include_mono)) |> tidyr::pivot_longer(paste0(calib_samples_mono), names_to ="calib", values_to = paste0("calib_mono.xlsx"))
  
for (name in names(paradise_reports_list)) {
  data <- paradise_reports_list[[name]] |> dplyr::filter(New_CAS_Name %in% c("(±)-α-Pinene", "β-Pinene", "Limonene", "p-Cymene (8CI)" ))|>
  dplyr::select(tidyselect::all_of(cols_to_include_mono)) |> tidyr::pivot_longer(paste0(calib_samples_mono), names_to ="calib", values_to = paste0(name))


  calib_btw_batch_mono <- dplyr::full_join(calib_btw_batch_mono, data, by = c("calib", "New_CAS_Name",  "Hit 1: CAS" ))
}
  calib_btw_batch_mono <- dplyr::full_join(calib_btw_batch_mono, calib_solo_mono, by = c("calib", "New_CAS_Name",  "Hit 1: CAS" ))

  
  return(  calib_btw_batch_mono)
}
