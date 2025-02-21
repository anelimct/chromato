
plot_in_out <- function(time, values_T_in, values_T_out, values_H_in, values_H_out, name, year){
  df <- data.frame(
  Time = as.POSIXct(unlist(time), format = "%Y-%m-%d %H:%M:%S", tz = "CET"),
  T_in = unlist(values_T_in),
  T_out = unlist(values_T_out),
  H_in = unlist(values_H_in),
  H_out = unlist(values_H_out))
  
  df <- df |> tidyr::pivot_longer( cols = c("T_in", "T_out", "H_in", "H_out" ), names_to = "variable", values_to = "value" )
  
  t_title <- rlang::englue("Temperature for {name}")
  h_title <- rlang::englue("Humidity for {name}")
  
  
  gg_t <- ggplot(df |> dplyr::filter(variable %in% c("T_in", "T_out")),
                    aes(x = Time, y = value, color = variable)) +
    geom_line(size = 1) +
    labs(
         x = "Temps",
         y = "Température (°C)",
         color = "variable") +
    ylim(18, 55)+
    theme_minimal()+
    ggtitle(t_title)
  
  gg_h <- ggplot(df |> dplyr::filter(variable %in% c("H_in", "H_out")),
                    aes(x = Time, y = value, color = variable)) +
    geom_line(size = 1) +
    labs(
      x = "Temps",
      y = "Humidité (%RH)",
      color = "variable") +
    theme_minimal()+
    ggtitle(h_title)
  
  gg_final <- cowplot::plot_grid(gg_t, gg_h)
  plot_path <- stringr::str_glue("figures/graphs_ibuttons/{year}")
  
  
  folder_created = {
    if (!dir.exists(plot_path)) {
      dir.create(plot_path, recursive = TRUE)
    }
    plot_path
  }
  
  
  plot_name <- stringr::str_glue("{tolower(name)}_ibuttons.png")
  
  complete_plot_path <- stringr::str_glue("figures/graphs_ibuttons/{tolower(name)}_ibuttons.png")
  
    ggsave( plot_name, gg_final, path = plot_path, width = 16, height = 4.95, units = "cm", bg = "white" )
  return(complete_plot_path)
}

save_plot_ibuttons <- function (data, year){
  
  # enlever les lignes pour lesquelles on n'a pas toutes les données 
  data <- data |> 
    dplyr::filter(
      lengths(Time_T_in) > 0 &
        lengths(values_T_in) > 0 &
        lengths(values_T_out) > 0 &
        lengths(values_H_in) > 0 &
        lengths(values_H_out) > 0
    )
  
  
  purrr::pmap(list( time = data$Time_T_in, values_T_in = data$values_T_in, values_T_out = data$values_T_out, values_H_in = data$values_H_in, values_H_out = data$values_H_out, name = data$ID , year = year), plot_in_out)
}


