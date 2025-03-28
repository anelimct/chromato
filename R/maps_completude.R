# 
# WOODIV_grid <- sf::st_read(paste0(here::here ( "data", "WOODIV", "SPATIAL", "WOODIV_grid"), "/WOODIV_grid.shp"))     # Remplacez par le chemin correct
# WOODIV_shape <- sf::st_read(paste0(here::here ( "data", "WOODIV", "SPATIAL", "WOODIV_shape"), "/WOODIV_shape.shp"))
# load(paste0(here::here ( "data", "WOODIV"), "/working_file.rdata"))




compute_completeness <- function(WOODIV_grid, working_file, summary_DB) {
  # Join the datasets to create WOODIV_data
  WOODIV_data <- WOODIV_grid |>
    dplyr::left_join(working_file, by = c("Idgrid" = "Idgrid")) |>
    dplyr::left_join(summary_DB, by = c("to_aggregate_with" = "spcode.agg"))
  
  # Calculate completeness
  completeness <- WOODIV_data |>
    dplyr::group_by(Idgrid) |>
    dplyr::summarise(
      total_species = dplyr::n_distinct(to_aggregate_with),
      bvocs_completeness = dplyr::n_distinct(to_aggregate_with[!is.na(isoprene)]) / total_species * 100
    ) |>
    dplyr::ungroup()
  
  return(completeness)
}


map_et_plot_completness <- function (completeness, WOODIV_shape){
  
  map <- ggplot() +
    geom_sf(data = WOODIV_shape, fill = "lightgrey", color = NA) +  # Add the background grid
    geom_sf(data = completeness, aes(fill = bvocs_completeness), size = 0.5) +  # Add the completeness data
    scale_fill_gradient(low = "white", high = "#03a219", name = "Completeness (%)") +
    theme_minimal() +
    labs(title = "Map of Data Completeness",
         fill = "Completeness") +
    theme(
      panel.grid.major = element_blank(),  # Remove major grid lines
      panel.grid.minor = element_blank(),  # Remove minor grid lines
      legend.position = "bottom"
    )
  
  hist <- ggplot(completeness, aes(x = bvocs_completeness)) +
    geom_histogram(binwidth = 10, fill = "#d20a21", color = "black", alpha = 0.7) +
    theme_minimal() +
    labs(title = "Distribution of Data Completeness",
         x = "Completeness (%)",
         y = "Frequency") +
    theme(
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank()
    )
  
  ggsave( "map_completness_bvocs.png",   map, path = paste0(here::here ("figures"),"map_completness_bvocs.png" ), width = 16, height =  16, units = "cm", bg = "white" )
  
  ggsave( "hist_completness_bvocs.png",   hist, path = paste0(here::here ("figures"),"hist_completness_bvocs.png" ), width = 16, height =  16, units = "cm", bg = "white" )
  
  return(completeness)
}

ranking_species <- function(working_file){
  data <- working_file |>  dplyr::rename( spcode.agg = to_aggregate_with ) |>  dplyr::mutate(Taxon = paste0(genus, " ", species)) |> dplyr::group_by(spcode.agg, Taxon) |> dplyr::summarise(nb_grid = dplyr::n_distinct(Idgrid)) |>  dplyr::ungroup()|> dplyr::mutate(rank = rank(-nb_grid))
}

