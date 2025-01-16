make_ibuttons_table <- function(paths_csv_files, type){
  # Keep only csv files on the type of variables retained
  paths_csv_files <- grep( paste0("*\\_", type, ".csv"), paths_csv_files, value= TRUE)
  
  #découper le path en autant d'élément que de / et gader le dernier qui est le nom du fichier
  name_file <- stringr::str_split(paths_csv_files, "/") |> 
  lapply( function(x){
    x[length(x)]
    })
  
  N_buttons <- stringr::str_sub(name_file, end = 2)
  # mettre des "-" entre les différent élément de la date qui était dans le nom du fichier ou tout était collé
  dates <- stringr::str_sub(name_file, start = 4, end = 11) |> 
    (\(x) stringr::str_c(
      stringr::str_sub(x, end = 4),
      stringr::str_sub(x, start = 5, end = 6),
      stringr::str_sub(x, start = 7, end = 8),
      sep = "-"
    ))()
    
  # read all the files use fill = T because in some files there is missing the extra value at this end and no "," so otherwise it stops
  list.DFs <- lapply(paths_csv_files, function(file) {
    data.table::fread(file, fill = T, skip = 19)
  })

  #coller le numéro du ibutton et la date qui proviennent du nom du fichier
  list.DFs <- mapply(function(df, n) {
    dplyr::mutate(df, N_button = n)
  }, list.DFs, N_buttons, SIMPLIFY = FALSE)

  list.DFs <- mapply(function(df, n) {
    dplyr::mutate(df, Date_file_name = n)
  }, list.DFs, dates, SIMPLIFY = FALSE)
  
  
  df <- do.call(rbind, list.DFs) |>
    dplyr::mutate( Date_Time = stringr::str_replace( `Date/Time`, "/24", "/2024")) |> #day ne peut jamais avoir un slash devant puisque en première position donc c'est OK
    dplyr::mutate( Date_Time = stringr::str_replace( Date_Time, "/23", "/2023")) |>
    dplyr::mutate(Date_Time = as.POSIXct(Date_Time, format = "%d/%m/%Y %H:%M:%OS"))
  
  # failed_file<- df |>  dplyr::mutate(Date_file_name = as.Date (Date_file_name, format = "%Y-%m-%d")) |> 
  #   dplyr::filter(Date_file_name != as.Date(Date_Time))
  #   if (nrow(failed_file) > 0) {
  #     warning("Certains nom de fichiers ne correspondent pas à la date interne")}
  
  ## Pour les fichiers temp tous les unique(fail[,c("Date_file_name", "N_button")]) c'est la date du premier jour qui est indiqué sauf 44_20210724 où c'est le deuxième
    return(df)
}

export_ibuttons_data <- function(data, ibutton_table, type){
  variable <- paste0("ibutton_", type)
  
  #nettoyer le data pour qu'il y est deux digits aux N..iButton.IN et N..iButton.OUT et s'assurer que le h de heure est en minuscule
  data <- data |>
    dplyr::mutate(N..iButton.IN = stringr::str_pad(N..iButton.IN, 2, side = "left", pad = "0")) |> 
    dplyr::mutate(N..iButton.OUT = stringr::str_pad(N..iButton.OUT, 2, side = "left", pad = "0")) |> 
    dplyr::mutate(Heure.début = stringr::str_to_lower(Heure.début)) |>
    dplyr::mutate(Heure.fin = stringr::str_to_lower(Heure.fin))
  
  
  #faire des intervalles pour les périodes de prélevement 
  data <- data |>
    dplyr::mutate(start_prelevement = stringr::str_c(Date, Heure.début, sep=" ")) |> 
    dplyr::mutate(end_prelevement = stringr::str_c(Date, Heure.fin , sep=" ")) |> 
    dplyr::mutate(start_prelevement = as.POSIXct(start_prelevement, format = "%d/%m/%Y %Hh%M")) |> 
    dplyr::mutate(end_prelevement = as.POSIXct(end_prelevement, format = "%d/%m/%Y %Hh%M")) |> 
    dplyr::mutate(interval = lubridate::interval(start_prelevement, end_prelevement)) 
  
  # the !! paste0 is to have an interactive column name in the mutate in order to add the suffix _T ou _H to colname 
  data <- data |> 
    dplyr::rowwise() |> 
    dplyr::mutate(
      !!paste0("values_", type, "_in") := list(
        ibutton_table |> 
          dplyr::filter(Date_Time %within% interval & N_button == N..iButton.IN) |> 
          dplyr::pull(Value)
      ),
      !!paste0("values_", type, "_out") := list(
        ibutton_table |> 
          dplyr::filter(Date_Time %within% interval & N_button == N..iButton.OUT) |> 
          dplyr::pull(Value)
      )
    ) |> 
    dplyr::ungroup()
  return(data)    
      
}
