# library(elevatr)
# library(terra)
# library(exactextractr)
# library(sf)
# 
# # Télécharger le MNT à résolution ~1 km (z=8)
# elev_raster <- get_elev_raster(
#   locations = WOODIV_wgs84,
#   z = 8,                     # Résolution ~1 km
#   src = "aws",                # Données SRTM AWS (pas de clé)
#   clip = "tile",         # Découper à l'emprise des grilles
#   override_size_check = TRUE  # Accepte la taille estimée (modérée)
# )
