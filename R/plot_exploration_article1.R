plot_pie_chart <- function(data, type_column = "type", title = "Distribution des types") {
  # Charger les packages nécessaires
  library(plotly)
  library(dplyr)
  
  # Préparer les données : compter les occurrences de chaque type
  data_summary <- data |> 
    dplyr::rename(type = !!rlang::sym(type_column)) |> 
    dplyr::count(type, name = "count") |> 
    dplyr::arrange(desc(count))
  
  # Mapping des labels
  label_mapping <- c(
    "mono" = "Monoterpenes",
    "iso" = "Isoprene", 
    "both" = "Both",
    "NE" = "Non-emitters"
  )
  
  # Mapping des couleurs
  color_mapping <- c(
    "iso" = "#70c230",
    "both" = "#da6484",
    "mono" = "#3764da", 
    "NE" = "#fce72e"
  )
  
  # Appliquer le mapping des labels
  data_summary <- data_summary |> 
    dplyr::mutate(
      display_label = ifelse(type %in% names(label_mapping), 
                             label_mapping[type], 
                             type)
    )
  
  # Créer le vecteur de couleurs ordonné selon data_summary
  ordered_colors <- sapply(data_summary$type, function(x) {
    if (x %in% names(color_mapping)) {
      return(color_mapping[x])
    } else {
      # Couleur par défaut pour les types non spécifiés
      return("#cccccc")
    }
  })
  
  # Créer le graphique plotly
  fig <- plotly::plot_ly(
    data = data_summary,
    labels = ~display_label,
    values = ~count,
    type = 'pie',
    textposition = 'inside',
    textinfo = 'label+percent',
    insidetextfont = list(color = '#FFFFFF', size = 12),
    hoverinfo = 'text',
    hovertext = ~paste('<b>', display_label, '</b><br>',
                       'Type original: ', type, '<br>',
                       'Count: ', count, '<br>',
                       'Percentage: ', round(count/sum(count)*100, 1), '%'),
    marker = list(
      colors = ordered_colors,
      line = list(color = '#FFFFFF', width = 2)
    ),
    showlegend = TRUE,
    hole = 0  # Pie chart normal (0 pour pie chart, >0 pour donut chart)
  )
  
  # Personnaliser la mise en page
  fig <- fig |>  plotly::layout(
    title = list(
      text = title,
      font = list(size = 20, family = "Arial", color = "#333333")
    ),
    xaxis = list(
      showgrid = FALSE, 
      zeroline = FALSE, 
      showticklabels = FALSE
    ),
    yaxis = list(
      showgrid = FALSE, 
      zeroline = FALSE, 
      showticklabels = FALSE
    ),
    legend = list(
      title = list(text = "Types"),
      orientation = "v",
      font = list(size = 12),
      x = 1.05,
      y = 0.5
    ),
    margin = list(l = 20, r = 120, t = 50, b = 20)
  )
  
  return(fig)
}







create_density_facet_plot <- function(data,
                                      x_var = "SLA",
                                      group_var = "type",
                                      plot_type = "density_mean",
                                      same_y_scale = TRUE,
                                      show_histogram = TRUE,
                                      trans = "none") {   # nouveau paramètre
  
  library(ggplot2)
  library(dplyr)
  
  # Vérification des variables
  if (!x_var %in% names(data)) stop(paste("Variable", x_var, "non trouvée"))
  if (!group_var %in% names(data)) stop(paste("Variable", group_var, "non trouvée"))
  
  # Transformation de la variable
  data <- data %>%
    mutate(x_trans = case_when(
      trans == "log"  ~ log(.data[[x_var]]),
      trans == "sqrt" ~ sqrt(.data[[x_var]]),
      TRUE            ~ .data[[x_var]]
    ))
  
  # Couleurs et labels
  type_colors <- c(
    "iso"  = "#70c230",
    "both" = "#da6484",
    "mono" = "#3764da",
    "NE"   = "#fce72e"
  )
  
  type_labels <- c(
    "iso"  = "Isoprene",
    "both" = "Both",
    "mono" = "Monoterpenes",
    "NE"   = "Non-emitters",
    "All"  = "All species"
  )
  
  # ===== Créer les facettes =====
  data_facet <- data %>%
    mutate(.facet = as.character(.data[[group_var]]),
           .fill  = as.character(.data[[group_var]]))
  
  data_all <- data %>%
    mutate(.facet = "All",
           .fill  = "All")
  
  data_plot <- bind_rows(data_facet, data_all)
  
  # ===== Calcul des moyennes =====
  if (plot_type %in% c("density_mean", "combined")) {
    mean_groups <- data %>%
      group_by(.data[[group_var]]) %>%
      summarise(mean_value = mean(x_trans, na.rm = TRUE), .groups = "drop") %>%
      mutate(.facet = as.character(.data[[group_var]]),
             .fill  = as.character(.data[[group_var]]))
    
    mean_all <- data %>%
      summarise(mean_value = mean(x_trans, na.rm = TRUE)) %>%
      mutate(.facet = "All",
             .fill  = "All")
    
    group_means <- bind_rows(mean_groups, mean_all)
  }
  
  # ===== Graphique de base =====
  p <- ggplot(data_plot, aes(x = x_trans, fill = .fill))
  
  # Histogramme
  if (show_histogram) {
    p <- p +
      geom_histogram(aes(y = after_stat(density)),
                     bins = 30,
                     color = "white",
                     alpha = 0.4)
  }
  
  # Densité
  p <- p + geom_density(alpha = 0.6, size = 0.8, color = "black")
  
  # Lignes des moyennes
  if (plot_type %in% c("density_mean", "combined")) {
    p <- p +
      geom_vline(data = group_means,
                 aes(xintercept = mean_value, color = .fill),
                 linetype = "dashed",
                 size = 1)
  }
  
  # Facettes
  p <- p +
    facet_grid(.facet ~ .,
               scales = if (same_y_scale) "fixed" else "free_y",
               labeller = as_labeller(type_labels))
  
  # Couleurs manuelles
  p <- p +
    scale_fill_manual(values = c(type_colors, "All" = "grey70"), guide = "none") +
    scale_color_manual(values = c(type_colors, "All" = "black"), guide = "none")
  
  # ===== Labels des axes avec transformation =====
  x_label <- switch(trans,
                    "log"  = paste0("log(", x_var, ")"),
                    "sqrt" = paste0("sqrt(", x_var, ")"),
                    x_var)
  
  p <- p +
    labs(
      title = paste("Distribution de", x_label),
      x = x_label,
      y = "Densité"
    ) +
    theme_minimal() +
    theme(
      strip.text = element_text(face = "bold"),
      strip.background = element_rect(fill = "grey95"),
      plot.title = element_text(hjust = 0.5, face = "bold")
    )
  
  return(p)
}




# Charger les packages nécessaires
library(corrplot)
library(ggplot2)
library(dplyr)



matrice_correlation_spearmann <- function(data, include_type = FALSE) {
  
  # Sélectionner les traits d'intérêt
  traits <- c("HeightMax", "LeafArea", "SLA", "StemSpecDens")
  
  # Inclure la variable "type" si demandé
  if (include_type) {
    traits <- c(traits, "type")
  }
  
  # Créer un sous-ensemble des données
  data_traits <- data[, traits]
  
  # Matrice de corrélation de Spearman pour TOUTES les espèces
  M <- cor(data_traits[, 1:4], 
           method = "spearman")
  
  
  
  col <- colorRampPalette(c("#BB4444", "#EE9988", "#FFFFFF", "#77AADD", "#4477AA"))
  
  
  p.mat <- cor.mtest(data_traits[, 1:4], method = "spearman")
  
  
  

  par(mfrow = c(1, 2))
  
  # Corrplot
  corrplot(M,
           method = "color",
           type = "upper",
           addCoef.col = "black",
           p.mat = p.mat$p,
           sig.level = 0.05,
           insig = "blank")
  
  title("Spearman correlations")
  
  # Scatterplot matrix
  pairs(data_traits[, 1:4],
        pch = 19,
        col = rgb(0, 0, 0, 0.5),
        main = "Scatterplot matrix", lower.panel = NULL, panel = panel.smooth)
  

}


matrice_correlation_spearmann_by_type <- function(data){
  # Vérifier les types disponibles
  types <- unique(data$type)
  print(paste("Types disponibles:", paste(types, collapse = ", ")))
  
  
  for (t in types) {
    # Sous-ensemble pour le type t
    subset_data <- data |> dplyr:: filter(type == t)
    
    
    if(nrow(subset_data) < 8){
      warning(paste("Type", t, ": pas assez de données → ignoré"))
      next
    }
    
    
    matrice_correlation_spearmann(subset_data)
    title(main = paste("Spearman – Type:", t))
    
  
}
}


#test <- merged_data_3T|>  subset( select = c(isoprene, monoterpenes, HeightMax, LeafArea, SLA, SeedMass, StemSpecDens)) |> normaliser_dataframe()
#clean_test <- na.omit(test[, c("isoprene", "monoterpenes", "HeightMax", "LeafArea", "SLA", "SeedMass", "StemSpecDens")])







plot_funspace_with_species_and_spearman <- function(
    pca_obj,  # Objet PCA (ex: pca.trait_usual_iso)
    group_vec,  # Vecteur de groupes (ex: type_isoprenoids$type)
    noms_a_afficher,  # Noms des espèces à afficher
    main_title,  # Titre du plot
    data_spearman,  # Données pour le test de Spearman
    PCs = c(1, 2)  # Composantes principales à afficher
) {
  
  # 1. Créer l'objet funspace
  trait_space <- funspace::funspace(pca_obj, PCs = PCs, group.vec = group_vec)
  
  # 2. Extraire les coordonnées des espèces
  species_coords <- trait_space$originalX$scores
  
  # 3. Plot funspace avec les espèces et leurs étiquettes
  plot(trait_space, quant.plot = TRUE, arrows = TRUE, arrows.length = 0.9)
  
  # Extraire les coordonnées des espèces à afficher
  coords_a_afficher <- species_coords[rownames(species_coords) %in% noms_a_afficher, c("Comp.1", "Comp.2")]
  
  # Tracer les points précis pour chaque espèce
  points(coords_a_afficher[, "Comp.1"], coords_a_afficher[, "Comp.2"],
         pch = 19, col = "black", cex = 1.2)
  
  # Appliquer un jitter global aux coordonnées des étiquettes
  coords_jittered <- jitter(as.matrix(coords_a_afficher), amount = 0.02)
  
  # Ajouter les étiquettes avec le jitter global
  text(coords_jittered[, 1], coords_jittered[, 2],
       labels = rownames(coords_jittered),
       cex = 0.8, col = "black", pos = 3)
  
  # Ajouter un titre
  title(main = main_title, cex.main = 1.2)
  
  # 4. Plot de type "groups"
  plot(x = trait_space, type = "groups", quant.plot = TRUE, globalContour = TRUE,
       pnt = TRUE, pnt.cex = 0.1, pnt.col = rgb(0.2, 0.8, 0.1, alpha = 0.2),
       axis.title.line = 1)
  
  # 5. Test de Spearman et corrélations
  matrice_correlation_spearmann(data_spearman)
}





