# 
# WOODIV_grid <- sf::st_read(paste0(here::here ( "data", "WOODIV", "SPATIAL", "WOODIV_grid"), "/WOODIV_grid.shp"))     # Remplacez par le chemin correct
# WOODIV_shape <- sf::st_read(paste0(here::here ( "data", "WOODIV", "SPATIAL", "WOODIV_shape"), "/WOODIV_shape.shp"))
# load(paste0(here::here ( "data", "WOODIV"), "/working_file.rdata"))

process_summary_data <- function(summary_data, working_file) {
  final_data <- summary_data |>
    dplyr::mutate(
      all_distinct_origins_isoprene = dplyr::coalesce(distinct_origins_isoprene, 0) + dplyr::coalesce(distinct_origins, 0),
      all_distinct_origins_monoterpenes = dplyr::coalesce(distinct_origins_monoterpenes, 0) + dplyr::coalesce(distinct_origins, 0),
      all_nb_trees_isoprene = dplyr::coalesce(nb_trees_isoprene, 0) + dplyr::coalesce(n_entries, 0),
      all_nb_trees_monoterpenes = dplyr::coalesce(nb_trees_monoterpenes, 0) + dplyr::coalesce(n_entries, 0),
      all_nb_entries_isoprene = dplyr::coalesce(nb_entries_isoprene, 0) + dplyr::coalesce(n_entries, 0),
      all_nb_entries_monoterpenes = dplyr::coalesce(nb_entries_monoterpenes, 0) + dplyr::coalesce(n_entries, 0), 
      relative_grid = (nb_grid / length(unique(working_file$Idgrid)))*100
    )
  return(final_data)
}


compute_completeness <- function(WOODIV_grid, working_file, summary_all, minimum_pop) {
  # Join the datasets to create WOODIV_data
  WOODIV_data <- WOODIV_grid |>
    merge(working_file, by = c("Idgrid" = "Idgrid")) |>
    dplyr::left_join(summary_all, by = c("to_aggregate_with" = "spcode.agg"))
  
  # Calculate completeness
  completeness <- WOODIV_data |>
    dplyr::group_by(Idgrid, geometry) |>
    dplyr::summarise(
      total_species = dplyr::n_distinct(to_aggregate_with),
      bvocs_completeness_litt = dplyr::n_distinct(to_aggregate_with[!is.na(distinct_origins_isoprene) & !is.na(distinct_origins_monoterpenes) & distinct_origins_isoprene >= as.numeric(minimum_pop) & distinct_origins_monoterpenes >= as.numeric(minimum_pop)], na.rm = TRUE) / total_species * 100,
      bvocs_completeness_all = dplyr::n_distinct(to_aggregate_with[!is.na(all_distinct_origins_isoprene) & !is.na(all_distinct_origins_monoterpenes) & all_distinct_origins_isoprene >= as.numeric(minimum_pop) & all_distinct_origins_monoterpenes >= as.numeric(minimum_pop)], na.rm = TRUE) / total_species * 100
    ) |>
    dplyr::ungroup()
  
  return(completeness)
}


map_et_plot_completness <- function (completeness, WOODIV_shape){
  
  map <- ggplot() +
    geom_sf(data = WOODIV_shape, fill = "lightgrey", color = NA) +  # Add the background grid
    geom_sf(data = completeness, aes(fill = bvocs_completeness_litt), size = 0.5) +  # Add the completeness data
    scale_fill_gradient(low = "white", high = "#03a219", name = "Completeness (%)") +
    theme_minimal() +
    labs(title = "Map of Data Completeness from Litterature",
         fill = "Completeness") +
    theme(
      panel.grid.major = element_blank(),  # Remove major grid lines
      panel.grid.minor = element_blank(),  # Remove minor grid lines
      legend.position = "bottom"
    )
  
  map_2 <- ggplot() +
    geom_sf(data = WOODIV_shape, fill = "lightgrey", color = NA) +  # Add the background grid
    geom_sf(data = completeness, aes(fill = bvocs_completeness_all), size = 0.5) +  # Add the completeness data
    scale_fill_gradient(low = "white", high = "#03a219", name = "Completeness (%)") +
    theme_minimal() +
    labs(title = "Map of Data Completeness after field campaigns 2023_2024",
         fill = "Completeness") +
    theme(
      panel.grid.major = element_blank(),  # Remove major grid lines
      panel.grid.minor = element_blank(),  # Remove minor grid lines
      legend.position = "bottom"
    )
  
  hist <- ggplot() +
    geom_histogram(data = completeness, aes(x = bvocs_completeness_litt),
                  binwidth = 5, fill = "#d20a21", color = "white", alpha = 0.3) +
    geom_histogram(data = completeness, aes(x = bvocs_completeness_all),
                   binwidth = 5, fill = "#d20a21", color = "white", alpha = 0.7)+
    theme_minimal() +
    labs(title = "Distribution of Data Completeness",
         x = "Completeness (%)",
         y = "Frequency") +
    theme(
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank()
    )
  
  ggsave( "map_completness_bvocs_litt.png",   map, path = paste0(here::here ("figures", "completness")), width = 16, height =  16, units = "cm", bg = "white", create.dir = TRUE )
  ggsave( "map_completness_bvocs_all.png",   map_2, path = paste0(here::here ("figures", "completness")), width = 16, height =  16, units = "cm", bg = "white", create.dir = TRUE )
  
  ggsave( "hist_completness_bvocs_evolution.png",   hist, path = paste0(here::here ("figures", "completness")), width = 16, height =  16, units = "cm", bg = "white", create.dir = TRUE )
  
  return(completeness)
}

ranking_species <- function(working_file){
  data <- working_file |>  dplyr::rename( spcode.agg = to_aggregate_with ) |>  dplyr::mutate(Taxon = paste0(genus, " ", species)) |> dplyr::group_by(spcode.agg, Taxon) |> dplyr::summarise(nb_grid = dplyr::n_distinct(Idgrid)) |>  dplyr::ungroup()|> dplyr::mutate(rank = rank(-nb_grid))
}


plot_hist_ranking <- function(data, column_to_plot, tronquer_min, main_lab) {
  # Ordonner les données par relative_grid
  data_ordered <- data |>  
    dplyr::distinct_at("spcode.agg", .keep_all = T) |> 
    dplyr::arrange(relative_grid) |> 
    dplyr::mutate(Taxon.x = factor(Taxon.x, levels = unique(Taxon.x))) |> 
    dplyr::mutate(!!column_to_plot := as.factor(.data[[column_to_plot]])) |> 
    dplyr::filter(relative_grid >= as.numeric(tronquer_min))
  
  ggplot(data_ordered, aes(x = relative_grid, y = Taxon.x, fill = .data[[column_to_plot]])) +
    geom_bar(stat = "identity") +
    scale_fill_manual(values = c("white", "#f0b4c1","#e48da1", "#d2205f", "#a10f3b","#73054b", "#500334"  )) +
    labs(x = "Présence (%)", y = "Species Name", fill = "Number of population sampled", title = main_lab) +
    theme_minimal()
  
  ggsave(paste0(column_to_plot, ".png"), plot = last_plot(), path = paste0(here::here ("figures", "completness", "effort_echantillonnage")) , width = 3048, height = 2095, create.dir = T, limitsize = FALSE, units = "px", bg = "white")
  
  return(data)
}

plot_tree_effort_ech <- function ( data, tree){
  data_2 <- data |> 
    dplyr::mutate(all_distinct_origins_isoprene = ifelse(is.na(all_distinct_origins_isoprene), 0, all_distinct_origins_isoprene)) |> 
    dplyr::mutate(all_distinct_origins_monoterpenes = ifelse(is.na(all_distinct_origins_monoterpenes), 0, all_distinct_origins_monoterpenes)) |> 
    dplyr::mutate(all_distinct_origins_isoprene = as.factor(all_distinct_origins_isoprene)) |> 
    dplyr::mutate(all_distinct_origins_monoterpenes = as.factor(all_distinct_origins_monoterpenes)) |> 
    dplyr::mutate(Taxon.x = stringr::str_replace(Taxon.x, " ", "_")) |> 
    dplyr::select(Taxon.x, everything())
  
  

  
  # Prune the tree to keep only the specified tips
  p_tree<- ggtree::ggtree(tree)
  
  p <- ggtree::`%<+%`( p_tree, data_2) 
  p_iso <- p + ggtree::geom_tippoint(aes(color= all_distinct_origins_isoprene)) +
    scale_color_manual(values = c("white", "#f0b4c1", "#e48da1", "#d2205f", "#a10f3b", "#73054b", "#500334"))
  
  p_mono <- p + ggtree::geom_tippoint(aes(color= all_distinct_origins_monoterpenes)) +
    scale_color_manual(values = c("white", "#f0b4c1", "#e48da1", "#d2205f", "#a10f3b", "#73054b", "#500334"))
  
  ggsave("chronogram_iso.png", plot = p_iso, path = paste0(here::here ("figures", "completness", "effort_echantillonnage")) , width = 3048, height = 2095, create.dir = T, limitsize = FALSE, units = "px", bg = "white")
  ggsave("chronogram_monoterpenes.png", plot = p_mono, path = paste0(here::here ("figures", "completness", "effort_echantillonnage")) , width = 3048, height = 2095, create.dir = T, limitsize = FALSE, units = "px", bg = "white")
  
  return(data)
}
