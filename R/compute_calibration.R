compare_calib_btw_reports <- function (calib_files, paradise_reports_list, paradise_reports_calib_list , type){
 
  #puisque les samples sont en .CDF et non D les renommer
  calib_files <- calib_files |> dplyr::mutate(new_name =stringr::str_replace(new_name, ".D", ".CDF") )
  
  #liste de tous les noms des calib mono 
  calib_samples_mono <- unique(calib_files[stringr::str_detect( calib_files$new_name, "mono"), 1])
  
  calib_samples_iso <- unique(calib_files[stringr::str_detect( calib_files$new_name, "iso"), 1])
  #Toutes les colonnes que je veux récupérer
  cols_to_include_mono <- c(calib_samples_mono,  "New_CAS_Name", "Hit 1: CAS" )
  cols_to_include_iso <- c(calib_samples_iso,  "New_CAS_Name", "Hit 1: CAS" )
  
  #initialisationd d'un dataframe
  calib_btw_batch <- paradise_reports_list[[1]] |> dplyr::select("New_CAS_Name", "Hit 1: CAS") |> tibble::add_column(calib = as.character(NA)) 
  
  calib_btw_batch_mono <- calib_btw_batch[0,]
  calib_btw_batch_iso <- calib_btw_batch[0,]
  
  #Pur les calib mono aire mesurées dans une paradise session aves que des chromato calib, 
  #on ne garde que les lignes qui correspondent au mélange de mono utilisés
  #garder que les colonnes qui nous interesse CaD , les colonnes calib, hit CAS, et Name CAS
  #pivoter pour avoir le nom de la calib dans une colonne et les valeurs dans une colonne qui porte le nom de la paradise session
  calib_solo_mono <- paradise_reports_calib_list[["calib_mono.xlsx"]] |> dplyr::filter(New_CAS_Name %in% c("(±)-α-Pinene", "β-Pinene", "D-Limonene", "p-Cymene" ))|>
    dplyr::select(tidyselect::all_of(cols_to_include_mono)) |> tidyr::pivot_longer(paste0(calib_samples_mono), names_to ="calib", values_to = paste0("calib_mono.xlsx"))
  
  calib_solo_iso <- paradise_reports_calib_list[["calib_iso.xlsx"]] |> dplyr::filter(New_CAS_Name %in% c("Isoprene" ))|>
    dplyr::select(tidyselect::all_of(cols_to_include_iso)) |> tidyr::pivot_longer(paste0(calib_samples_iso), names_to ="calib", values_to = paste0("calib_iso.xlsx"))
  
if (type == "mono") {
  
  for (name in names(paradise_reports_list)) {
    data <- paradise_reports_list[[name]] |> dplyr::filter(New_CAS_Name %in% c("(±)-α-Pinene", "β-Pinene", "D-Limonene", "p-Cymene" ))|>
      dplyr::select(tidyselect::all_of(cols_to_include_mono)) |> tidyr::pivot_longer(paste0(calib_samples_mono), names_to ="calib", values_to = paste0(name))
    
    
    calib_btw_batch_mono <- dplyr::full_join(calib_btw_batch_mono, data, by = c("calib", "New_CAS_Name",  "Hit 1: CAS" ))
  }
  calib_btw_batch_mono <- dplyr::full_join(calib_btw_batch_mono, calib_solo_mono, by = c("calib", "New_CAS_Name",  "Hit 1: CAS" ))
  result <- calib_btw_batch_mono
} else {
  
  
  for (name in names(paradise_reports_list)) {
    data <- paradise_reports_list[[name]] |> dplyr::filter(New_CAS_Name %in% c("Isoprene" ))|>
      dplyr::select(tidyselect::all_of(cols_to_include_iso)) |> tidyr::pivot_longer(paste0(calib_samples_iso), names_to ="calib", values_to = paste0(name))
    
    
    calib_btw_batch_iso <- dplyr::full_join(calib_btw_batch_iso, data, by = c("calib", "New_CAS_Name",  "Hit 1: CAS" ))
  }
  calib_btw_batch_iso <- dplyr::full_join(calib_btw_batch_iso, calib_solo_iso, by = c("calib", "New_CAS_Name",  "Hit 1: CAS" ))
  result <-   calib_btw_batch_iso
}
  
  
  return( result)
}




plot_calib_btw_session <- function(table_calib_btw_session, calib_file, type = c("iso", "mono")) {
  type <- match.arg(type)
  
  calib_file <- calib_file |>
    dplyr::select(new_name, µg) |>
    dplyr::mutate(new_name = stringr::str_replace(new_name, ".D", ".CDF")) |>
    dplyr::rename(calib = new_name)
  
  cols_to_plot <- colnames(table_calib_btw_session)
  cols_to_plot <- cols_to_plot[-c(1, 2, 3)]
  
  grouped_data <- dplyr::left_join(table_calib_btw_session, calib_file) |>
    dplyr::group_by(New_CAS_Name)
  
  # Loop through each group and create a plot and summary table
  for (name in unique(grouped_data$New_CAS_Name)) {
    # Filter data for the current group
    compound_data <- grouped_data |>
      dplyr::filter(New_CAS_Name == name)
    
    # Melt the data to long format for ggplot2
    long_data <- tidyr::pivot_longer(compound_data, cols = all_of(cols_to_plot), names_to = "session", values_to = "value")
    
    # Calculate model data
    model_data <- long_data |>
      dplyr::group_by(session) |>
      dplyr::summarise(
        slope = coef(lm(value ~ 0 + µg))[1],
        R2 = summary(lm(value ~ 0 + µg))$r.squared,
        p_value = summary(lm(value ~ 0 + µg))$coefficients[1, 4]
      )
    
    # Create a summary table
    summary_table <- model_data |>
      dplyr::select(session, slope, p_value, R2)
    
    # Create the plot
    p <- ggplot(long_data, aes(x = µg, y = value, color = session)) +
      geom_point() +  # Scatter plot
      theme_minimal() +
      labs(title = paste("Plot for", name), x = "µg", y = "Area") +
      geom_abline(
        data = model_data,
        aes(slope = slope, intercept = 0, color = session),
        linetype = "dashed"
      )
    
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
    compound_data <- grouped_data |>
      dplyr::filter(New_CAS_Name == name)
    
    for (i in seq_along(cols_to_plot)) {
      # Melt the data to long format for ggplot2
      long_data <- tidyr::pivot_longer(compound_data, cols = all_of(cols_to_plot[i]), names_to = "variable", values_to = "value")
      
      # Calculate the linear model
      model_data <- lm(value ~ 0 + µg, data = long_data)
      slope <- model_data$coefficients[1]
      
      p <- ggplot(long_data, aes(x = µg, y = value)) +
        geom_point() +
        theme_light() +
        labs(title = paste("Plot for", name, "batch", cols_to_plot[i]), x = "µg", y = "Area") +
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
  
  return(table_calib_btw_session)
}



