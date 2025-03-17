# 
# WOODIV_grid <- sf::st_read(paste0(here::here ( "data", "WOODIV", "SPATIAL", "WOODIV_grid"), "/WOODIV_grid.shp"))     # Remplacez par le chemin correct
# WOODIV_shape <- sf::st_read(paste0(here::here ( "data", "WOODIV", "SPATIAL", "WOODIV_shape"), "/WOODIV_shape.shp"))
# 
# 
# WOODIV_data <- WOODIV_grid %>%
#   left_join(working.file, by = c("Idgrid" = "Idgrid"))
# 
# 
# ggplot() +
#   geom_sf(data = WOODIV_shape, fill = "lightgrey", color = NA) +  # Ajouter la grille en arrière-plan
#   geom_sf(data = WOODIV_grid, size = 0.5) +   # Ajouter les occurrences
#   theme_minimal() +
#   labs(title = "Carte des Occurrences",
#        color = "Genus") +
#   theme(legend.position = "bottom")+
#   theme(
#     panel.grid.major = element_blank(),  # Supprimer les grilles principales
#     panel.grid.minor = element_blank(),  # Supprimer les grilles secondaires
#     legend.position = "bottom"
#   )

