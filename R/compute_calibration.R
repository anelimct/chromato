#' Title
#'
#' @param calib_quanti 
#' @param paradise_reports_list 
#' @param paradise_reports_calib_list 
#' @param type 
#'
#' @returns un tableau avec les composés, la calib ou on les a trouvé et l'aire associée a cette calib pour ce composé
#' @export
#'
#' @examples
compare_calib_btw_reports <- function (calib_quanti, paradise_reports_list, paradise_reports_calib_list , type){
  
  #puisque les samples sont en .CDF et non D les renommer
  calib_files <- calib_quanti |> dplyr::mutate(new_name =stringr::str_replace(new_name, ".D", ".CDF") )
  
  #liste de tous les noms des calib mono 
  calib_samples_mono <- unique(calib_files[stringr::str_detect( calib_files$new_name, "mono"), 1])
  
  calib_samples_iso <- unique(calib_files[stringr::str_detect( calib_files$new_name, "iso"), 1])
  
  #Toutes les colonnes que je veux récupérer
  cols_to_include_mono <- c(calib_samples_mono,  "New_CAS_Name" )
  cols_to_include_iso <- c(calib_samples_iso,  "New_CAS_Name" )
  
  
  calib_solo_iso <- paradise_reports_calib_list[["calib_iso.xlsx"]] |> dplyr::filter(New_CAS_Name %in% c("Isoprene" ))|>
    dplyr::select(tidyselect::all_of(cols_to_include_iso)) |> tidyr::pivot_longer(paste0(calib_samples_iso), names_to ="calib", values_to = paste0("calib_iso.xlsx"))
  
  
  #initialisationd d'un dataframe
  calib_btw_batch <- paradise_reports_list[[1]] |> dplyr::select("New_CAS_Name") |> tibble::add_column(calib = as.character(NA)) 
  
  calib_btw_batch_iso <- calib_btw_batch[0,]
  calib_btw_batch_mono <- calib_btw_batch[0,]
  
  #Pur les calib mono aire mesurées dans une paradise session aves que des chromato calib, 
  #on ne garde que les lignes qui correspondent au mélange de mono utilisés
  #garder que les colonnes qui nous interesse CaD , les colonnes calib, hit CAS, et Name CAS
  #pivoter pour avoir le nom de la calib dans une colonne et les valeurs dans une colonne qui porte le nom de la paradise session
  
  if (type == "mono") {
    
    for (name in names(paradise_reports_list)) {
      data <- paradise_reports_list[[name]] |> dplyr::select(tidyselect::any_of(cols_to_include_mono)) 
      
      calibartions <- intersect(calib_samples_mono, colnames(data))
      compounds <-   calib_files$compound[calib_files$new_name == calibartions[1]]
      
      data <- data |>  dplyr::filter(New_CAS_Name %in% compounds)|> tidyr::pivot_longer( cols = tidyselect::contains("mono"), names_to = "calib", values_to = name)
      
      
      calib_btw_batch_mono <- dplyr::full_join(calib_btw_batch_mono, data, by = c("calib", "New_CAS_Name" ))
    }
    result <- calib_btw_batch_mono
  } else {
    
    
    for (name in names(paradise_reports_list)) {
      data <- paradise_reports_list[[name]] |> dplyr::filter(New_CAS_Name %in% c("Isoprene" ))|>
        dplyr::select(tidyselect::any_of(cols_to_include_iso)) |>  tidyr::pivot_longer( cols = tidyselect::contains("iso"), names_to = "calib", values_to = name)
      
      
      calib_btw_batch_iso <- dplyr::full_join(calib_btw_batch_iso, data, by = c("calib", "New_CAS_Name" ))
    }
    calib_btw_batch_iso <- dplyr::full_join(calib_btw_batch_iso, calib_solo_iso, by = c("calib", "New_CAS_Name"))
    result <-   calib_btw_batch_iso
  }
  
  
  return( result)
}




plot_calib_btw_session <- function(table_calib_btw_session, calib_quanti, type = c("iso", "mono"), library_CAS_RI) {
  type <- match.arg(type)
  
  calib_file <- calib_quanti |>
    dplyr::select(new_name, compound, CAS, ng) |>
    dplyr::mutate(new_name = stringr::str_replace(new_name, ".D", ".CDF")) |>
    dplyr::rename(calib = new_name)
  
  cols_to_plot <- colnames(table_calib_btw_session)
  cols_to_plot <- cols_to_plot[-c(1, 2)]
  
  summary_table_all <- data.frame(
    session = character(),
    slope = numeric(),
    p_value = numeric(),
    R2 = numeric(),
    New_CAS_Name = character()
  )
  
  grouped_data <- dplyr::left_join(table_calib_btw_session, calib_file, by = c("calib" = "calib", "New_CAS_Name" = "compound")) |>
    dplyr::group_by(New_CAS_Name)
  
  # Loop through each group and create a plot and summary table
  for (name in unique(grouped_data$New_CAS_Name)) {
    # Filter data for the current group
    compound_data <- grouped_data |>
      dplyr::filter(New_CAS_Name == name)
    
    
    # Melt the data to long format for ggplot2
    long_data <- tidyr::pivot_longer(compound_data, cols = all_of(cols_to_plot), names_to = "session", values_to = "value") |> dplyr::filter(!is.na(value)) |> dplyr::distinct()
    
    # Calculate model data
    model_data <- long_data |>
      dplyr::group_by(session) |>
      dplyr::summarise(
        slope = coef(lm(value ~ 0+ng))[1],
        R2 = summary(lm(value ~ 0+ ng))$r.squared,
        p_value = summary(lm(value ~ 0+ng))$coefficients[1, 4]
      )
    
    # Create a summary table
    summary_table <- model_data |>
      dplyr::select(session, slope, p_value, R2)
    
    summary_table_longer <- summary_table |> dplyr::mutate(New_CAS_Name = name)
    
    summary_table_all <- rbind(summary_table_all , summary_table_longer ) 
    
    # Create the plot
    p <- ggplot(long_data, aes(x = ng, y = value, color = session)) +
      geom_point() +  # Scatter plot
      theme_minimal() +
      labs(title = paste("Plot for", name), x = "ng", y = "Area") +
      geom_abline(
        data = model_data,
        aes(slope = slope, intercept = 0, color = session),
        linetype = "dashed"
      )+ geom_hline(yintercept = 193500, color = "black") + 
      geom_hline(yintercept = 645000, color = "red")
    
    # Convert the summary table to a tableGrob for display
    table_grob <- gridExtra::tableGrob(summary_table, rows = NULL)
    
    # Arrange the plot and table side by side
    p_final <- gridExtra::grid.arrange(p, table_grob, ncol = 2)
    
    ggsave(plot = p_final, filename = paste0("plot_", name, "_all_session_", type, ".png"), path = here::here("figures", "graph_calib", paste0(type), "_all_sessions"), create.dir = TRUE, width = 35, height = 15, units = "cm")
  }

  grouped_data_2 <- grouped_data |>
    dplyr::ungroup() |>
    dplyr::group_by(New_CAS_Name, calib)

  for (name in unique(grouped_data$New_CAS_Name)) {
    # Filter data for the current group
    compound_data <- grouped_data_2 |>
      dplyr::filter(New_CAS_Name == name)

    for (i in seq_along(cols_to_plot)) {
      # Melt the data to long format for ggplot2
      long_data <- tidyr::pivot_longer(compound_data, cols = all_of(cols_to_plot[i]), names_to = "variable", values_to = "value")|> dplyr::filter(!is.na(value))
       if (length(long_data$value) == 0){
         
       } else{
         
         # Calculate the linear model
         model_data <- lm(value ~ 0 + ng, data = long_data)
         slope <- model_data$coefficients[1]
         
         p <- ggplot(long_data, aes(x = ng, y = value)) +
           geom_point() +
           theme_light() +
           labs(title = paste("Plot for", name, "batch", cols_to_plot[i]), x = "ng", y = "Area") +
           geom_abline(slope = slope, intercept = 0, linetype = "dashed")
         
         p2 <- ggplot() +
           geom_point(aes(x = model_data$fitted.values, y = model_data$residuals)) +
           labs(title = "Residuals vs Fitted", x = "Fitted values", y = "Residuals") +
           theme_light()
         
         p3 <- ggplot() +
           geom_qq(aes(sample = model_data$residuals)) +
           geom_qq_line(aes(sample = model_data$residuals)) +
           labs(title = "Normal Q-Q", x = "Theoretical Quantiles", y = "Sample Quantiles") +
           theme_light()
         
         bottom_plot <- cowplot::plot_grid(p2, p3, ncol = 2, labels = c('B', 'C'))
         # Combine the three plots into one window
         p_final <- cowplot::plot_grid(p, bottom_plot, ncol = 1, labels = c('A', ''))
         print(p_final)
         
         ggsave(plot = p_final, filename = paste0("plot_", name, "_", cols_to_plot[i], ".png"), path = here::here("figures", "graph_calib", paste0(type), "_by_session", paste0(name)), create.dir = TRUE, width = 15, height = 10, units = "cm")
         
       }
      
    }
  }


  summary_table_all<- dplyr::left_join(summary_table_all, library_CAS_RI) |>  
    dplyr::select("session", "slope", "p_value", "R2", "New_CAS_Name" , "New_CAS_Number" )
  return(summary_table_all)
}



