library(leaflet)
library(htmlwidgets)
library(here)

carte_notion_embed <- function(data, taxon_selectionne, nom_fichier = NULL) {
  
  if (is.null(nom_fichier)) {
    nom_fichier <- paste0("carte_", gsub(" ", "_", taxon_selectionne), ".html")
  }
  
  data_filtre <- data %>% 
    filter(Taxon == taxon_selectionne)
  
  if (nrow(data_filtre) == 0) {
    stop("Aucun prélèvement trouvé pour : ", taxon_selectionne)
  }
  
  # Créer carte avec aperçu optimisé
  carte <- leaflet(data_filtre, 
                   options = leafletOptions(
                     zoomControl = TRUE,
                     attributionControl = TRUE
                   )) %>%
    addTiles() %>%
    addMarkers(
      lng = ~Longitude..WGS84.,
      lat = ~Latitude..WGS84.,
      popup = ~paste(
        "<b>", Taxon, "</b><br>",
        "Date: ", Date, "<br>",
        "ID: ", ID, "<br>",
        "Lat: ", round(Latitude..WGS84., 4), "<br>",
        "Lon: ", round(Longitude..WGS84., 4)
      ),
      clusterOptions = markerClusterOptions()
    ) %>%
    setView(
      lng = mean(data_filtre$Longitude..WGS84.),
      lat = mean(data_filtre$Latitude..WGS84.),
      zoom = 8
    ) %>%
    addControl(
      paste("<div style='background: white; padding: 10px; border-radius: 5px;'>",
            "<h4>", taxon_selectionne, "</h4>",
            "<p>", nrow(data_filtre), " site(s) de prélèvement</p>",
            "</div>"),
      position = "topright"
    )
  
  # Sauvegarder avec métadonnées pour aperçu
  saveWidget(carte, 
             file = nom_fichier, 
             selfcontained = TRUE,
             title = paste("Carte", taxon_selectionne))
  
  message("Fichier créé : ", nom_fichier)
  message("\nPour Notion :")
  message("1. Uploader le fichier sur GitHub Pages, Netlify ou Google Drive")
  message("2. Copier le lien public")
  message("3. Utiliser le bloc 'Embed' dans Notion")
  
  return(carte)
}
