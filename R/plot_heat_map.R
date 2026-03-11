rename_mean_columns <- function(data, sp_screening) {
    # Créer un vecteur de correspondance spcode -> nom scientifique avec espaces
    name_mapping <- sp_screening |>
      dplyr::select(spcode, full_scientific_name) |>
      dplyr::distinct() |>
      dplyr::mutate(
        # Remplacer _ par espace
        clean_name = gsub("_", " ", full_scientific_name),
        # Nom de colonne original attendu (mean_ + spcode en minuscule)
        old_col_name = paste0("mean_", tolower(spcode))
      )

    # Créer un vecteur nommé pour le renommage
    rename_vector <- setNames(name_mapping$clean_name, name_mapping$old_col_name)

    # Renommer les colonnes
    for (i in seq_along(colnames(data))) {
      col_name <- colnames(data)[i]
      if (col_name %in% names(rename_vector)) {
        colnames(data)[i] <- rename_vector[[col_name]]
      }
    }

    return(data)
  }




my_compound_order <- c(
  "α-Terpinene",
  "α-Terpineol",
  "γ-Terpinene", 
  "Limonene",
  "Linalool",
  "β-Phellandrene",
  "2,3,6-Trimethyl-1,5-heptadiene (ACI)",
  "Cosmene",
  "Isoterpinolene",
 "Neoalloocimene",
 "Thuja-2,4(10)-diene",
 "alpha-Phellandrene", 
 "p-Cymenene",
 "p-Mentha-1,5,8-triene",
 "unknown monoterpene 1",
 "3-Carene",
 "β-Pinene",
 "α-Pinene",
 "Camphene",
 "β-Ocimene",
 "Isoprene",  
 "Terpinolene",
 "p-Cymene",
 "Sabinene",
 "1,8-Cineole",
 "Myrcene",
 "Camphor")




# Reverse the order
my_compound_order <- rev(my_compound_order)


my_species_order <- c(
  # Deciduous angiosperms
  "Betula pendula", "Betula pubescens", "Fraxinus angustifolia", "Fraxinus excelsior", "Ostrya carpinifolia",
  "Sorbus aucuparia", "Sorbus mougeotii", "Alnus viridis", "Salix eleagnos", "Salix pentandra", "Vitex agnus-castus",
  "Acer campestre", "Acer opalus", "Acer pseudoplatanus", "Alnus incana", "Corylus avellana", "Cotinus coggygria",
  "Ficus carica", "Frangula alnus", "Malus sylvestris", "Populus tremula",
  "Prunus avium", "Prunus brigantina", "Prunus mahaleb", "Prunus padus", "Pyrus pyraster", "Pyrus spinosa",
   "Rhamnus cathartica", "Salix alba", "Salix purpurea", "Salix cinerea", "Sambucus racemosa",
  "Sambucus nigra", "Sorbus aria", "Sorbus domestica", "Sorbus torminalis", "Tamarix africana", "Tamarix gallica",

  # Evergreen broadleaves (angiosperms)
  "Quercus crenata", "Ilex aquifolium", "Phillyrea angustifolia", "Rhamnus alaternus",

  # Evergreen conifers (gymnosperms)
  "Taxus baccata",
  "Juniperus thurifera",
  "Juniperus communis",
  "Pinus uncinata"
)


create_bubble_heatmap <- function(table, log_relative_prop = FALSE, 
                                  title = "",
                                  bubble_sizes = c(5, 50, 500),
                                  bubble_labels = NULL,
                                  compound_order = NULL,
                                  species_order = NULL) {
  
  # Transformer les données
  table_long <- table %>%
    tidyr::pivot_longer(
      cols = -compound,
      names_to = "species",
      values_to = "emission_rate"
    ) %>%
    dplyr::mutate(
      emission_rate_numeric = as.numeric(emission_rate),
      emission_rate_numeric = ifelse(is.na(emission_rate_numeric), 0, emission_rate_numeric)
    )
  
  # Calculer proportions
  table_long <- table_long %>%
    dplyr::group_by(species) %>%
    dplyr::mutate(
      total_species = sum(emission_rate_numeric, na.rm = TRUE),
      relative_proportion = ifelse(total_species > 0,
                                   (emission_rate_numeric / total_species) * 100,
                                   0)
    ) %>%
    dplyr::ungroup()
  
  # Filtrer
  table_filtered <- table_long %>%
    dplyr::filter(emission_rate_numeric > 0)
  
  # Appliquer l'ordre des composés si spécifié
  if(!is.null(compound_order)) {
    # Vérifier que tous les composés dans l'ordre existent dans les données
    existing_compounds <- unique(table_filtered$compound)
    compound_order <- compound_order[compound_order %in% existing_compounds]
    
    # Ajouter les composés manquants (s'ils ne sont pas dans l'ordre spécifié)
    missing_compounds <- setdiff(existing_compounds, compound_order)
    if(length(missing_compounds) > 0) {
      compound_order <- c(compound_order, missing_compounds)
    }
    
    # Convertir en facteur avec l'ordre spécifié
    table_filtered$compound <- factor(table_filtered$compound, levels = compound_order)
  } else {
    # Ordre par défaut : alphabétique
    table_filtered$compound <- factor(table_filtered$compound)
  }
  
  # Appliquer l'ordre des espèces si spécifié
  if(!is.null(species_order)) {
    # Vérifier que toutes les espèces dans l'ordre existent dans les données
    existing_species <- unique(table_filtered$species)
    species_order <- species_order[species_order %in% existing_species]
    
    # Ajouter les espèces manquantes
    missing_species <- setdiff(existing_species, species_order)
    if(length(missing_species) > 0) {
      species_order <- c(species_order, missing_species)
    }
    
    # Convertir en facteur avec l'ordre spécifié
    table_filtered$species <- factor(table_filtered$species, levels = species_order)
  } else {
    # Ordre par défaut : alphabétique
    table_filtered$species <- factor(table_filtered$species)
  }
  
  # Préparer les labels des bulles
  if(is.null(bubble_labels)) {
    bubble_labels <- paste0(bubble_sizes, " µg·g⁻¹·h⁻¹")
  }
  
  # Convertir les tailles de référence en échelle log
  reference_sizes <- log10(bubble_sizes + 1)
  
  if(log_relative_prop) {
    # Appliquer log aux proportions
    table_filtered <- table_filtered %>%
      dplyr::mutate(
        log_color = log10(relative_proportion + 1)
      )
    
    # Préparer échelle couleurs
    min_val <- min(table_filtered$log_color, na.rm = TRUE)
    max_val <- max(table_filtered$log_color, na.rm = TRUE)
    log_100 <- log10(100 + 1)
    if(max_val < log_100) max_val <- log_100
    
    color_breaks <- log10(c(1, 10, 100) + 1)
    color_labels <- c("1%", "10%", "100%")
    
    p <- ggplot2::ggplot(
      table_filtered,
      ggplot2::aes(
        x = species,
        y = compound,
        size = log10(emission_rate_numeric + 1),
        color = log_color
      )
    )
  } else {
    p <- ggplot2::ggplot(
      table_filtered,
      ggplot2::aes(
        x = species,
        y = compound,
        size = log10(emission_rate_numeric + 1),
        color = relative_proportion
      )
    )
  }
  
  p <- p +
    ggplot2::geom_point(alpha = 0.8) +
    
    # Échelle de taille personnalisée
    ggplot2::scale_size_continuous(
      name = expression("Emission factor" ~ (µg~.g^{-1}~.h^{-1})),
      range = c(2, 10),
      breaks = reference_sizes,
      labels = bubble_labels,
      guide = ggplot2::guide_legend(
        override.aes = list(
          color = "gray50",
          alpha = 0.8
        ),
        nrow = length(bubble_sizes),
        byrow = TRUE,
        label.position = "right",
        label.hjust = 0,
        keywidth = ggplot2::unit(1.5, "cm"),
        keyheight = ggplot2::unit(0.8, "cm")
      )
    )
  
  # Ajouter échelle de couleurs
  if(log_relative_prop) {
    p <- p + ggplot2::scale_color_viridis_c(
      name = "Log-transformed relative abundance (%)",
      option = "inferno",
      direction = -1,
      limits = c(min_val, max_val),
      breaks = color_breaks,
      labels = color_labels,
      guide = ggplot2::guide_colorbar(
        barwidth = 0.8,
        barheight = 12,
        ticks.colour = "white",
        frame.colour = "black"
      )
    )
  } else {
    p <- p + ggplot2::scale_color_viridis_c(
      name = "Relative proportion (%)",
      option = "inferno",
      direction = -1,
      limits = c(0, 100),
      breaks = c(0, 25, 50, 75, 100),
      guide = ggplot2::guide_colorbar(
        barwidth = 0.8,
        barheight = 12
      )
    )
  }
  
  # Thème final
  p <- p +
    ggplot2::theme_minimal() +
    ggplot2::theme(
      axis.text.x = ggplot2::element_text(angle = 45, hjust = 1, size = 8),
      axis.text.y = ggplot2::element_text(size = 7),
      panel.grid.major = ggplot2::element_line(color = "gray90"),
      legend.position = "right",
      legend.box = "vertical",
      legend.box.just = "left",
      plot.title = ggplot2::element_text(hjust = 0.5, size = 14, face = "bold"),
      legend.title = ggplot2::element_text(size = 10),
      legend.text = ggplot2::element_text(size = 8),
      legend.margin = ggplot2::margin(l = 10, r = 5),
      legend.spacing = ggplot2::unit(10, "pt")
    ) +
    ggplot2::labs(
      x = "Species",
      y = "Monoterpene compounds",
      title = title
    )
  
  return(p)
}
 library(dplyr)


create_bubble_heatmap_2 <- function(data,
                                    log_relative_prop = FALSE,
                                    title = "",
                                    bubble_sizes = c(1,10,50),
                                    bubble_labels = NULL,
                                    compound_order = NULL,
                                    species_order = NULL,
                                    isoprene_compound = "isoprene",
                                    bubble_sizes_isoprene = c(1,10,500),
                                    bubble_labels_isoprene = NULL,
                                    legend_position = "right",
                                    compound_groups = NULL,
                                    group_labels = NULL,
                                    bracket_left_offset = 0.35,
                                    bracket_tick_length = 0.12,
                                    bracket_text_offset = 0.25,
                                    save_plot = FALSE,
                                    file_name = "bubble_heatmap.png",
                                    dpi = 300,
                                    width = 210,
                                    height = 160) {
  
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(ggnewscale)
  
  #-----------------------------
  # COMPOUNDS THAT NEED *
  #-----------------------------
  
  starred_compounds <- c(
    "2,3,6-Trimethyl-1,5-heptadiene (ACI)",
    "Cosmene",
    "Isoterpinolene",
    "Neoalloocimene",
    "Thuja-2,4(10)-diene",
    "alpha-Phellandrene",
    "p-Cymenene",
    "p-Mentha-1,5,8-triene",
    "unknown monoterpene 1"
  )
  
  #-----------------------------
  # DATA PREPARATION
  #-----------------------------
  
  data_long <- data %>%
    pivot_longer(
      cols = -compound,
      names_to = "species",
      values_to = "emission_rate"
    ) %>%
    mutate(
      emission_rate_numeric = as.numeric(emission_rate),
      emission_rate_numeric = ifelse(is.na(emission_rate_numeric),0,emission_rate_numeric)
    )
  
  # Add "*" to selected compounds
  data_long <- data_long %>%
    mutate(
      compound = ifelse(compound %in% starred_compounds,
                        paste0(compound,"*"),
                        compound)
    )
  
  data_long <- data_long %>%
    group_by(species) %>%
    mutate(
      total_species = sum(emission_rate_numeric,na.rm=TRUE),
      relative_proportion = ifelse(total_species>0,
                                   emission_rate_numeric/total_species*100,
                                   0)
    ) %>%
    ungroup()
  
  data_filtered <- data_long %>%
    filter(emission_rate_numeric > 0)
  
  #-----------------------------
  # COMPOUND ORDER
  #-----------------------------
  
  if(!is.null(compound_order)){
    
    compound_order <- ifelse(compound_order %in% starred_compounds,
                             paste0(compound_order,"*"),
                             compound_order)
    
    existing_compounds <- unique(data_filtered$compound)
    
    compound_order <- compound_order[compound_order %in% existing_compounds]
    
    missing_compounds <- setdiff(existing_compounds,compound_order)
    
    compound_order <- c(compound_order,missing_compounds)
    
    data_filtered$compound <- factor(data_filtered$compound,
                                     levels=compound_order)
    
  } else {
    
    data_filtered$compound <- factor(data_filtered$compound)
    
  }
  
  #-----------------------------
  # SPECIES ORDER
  #-----------------------------
  
  if(!is.null(species_order)){
    
    existing_species <- unique(data_filtered$species)
    
    species_order <- species_order[species_order %in% existing_species]
    
    missing_species <- setdiff(existing_species,species_order)
    
    species_order <- c(species_order,missing_species)
    
    data_filtered$species <- factor(data_filtered$species,
                                    levels=species_order)
    
  } else {
    
    data_filtered$species <- factor(data_filtered$species)
    
  }
  
  # numeric coordinates
  
  data_filtered <- data_filtered %>%
    mutate(
      x_num = as.numeric(species),
      y_num = as.numeric(compound)
    )
  
  isoprene_data <- data_filtered %>% filter(compound %in% isoprene_compound)
  other_data <- data_filtered %>% filter(!(compound %in% isoprene_compound))
  
  #-----------------------------
  # BUBBLE SIZE LABELS
  #-----------------------------
  
  if(is.null(bubble_labels))
    bubble_labels <- paste0(bubble_sizes," µg·g⁻¹·h⁻¹")
  
  if(is.null(bubble_labels_isoprene))
    bubble_labels_isoprene <- paste0(bubble_sizes_isoprene," µg·g⁻¹·h⁻¹")
  
  reference_sizes <- log10(bubble_sizes+1)
  reference_sizes_isoprene <- log10(bubble_sizes_isoprene+1)
  
  #-----------------------------
  # COLOR TRANSFORMATION
  #-----------------------------
  
  if(log_relative_prop){
    
    data_filtered <- data_filtered %>%
      mutate(log_color = log10(relative_proportion+1))
    
    isoprene_data <- isoprene_data %>%
      mutate(log_color = log10(relative_proportion+1))
    
    other_data <- other_data %>%
      mutate(log_color = log10(relative_proportion+1))
    
    color_breaks <- log10(c(1,10,100)+1)
    color_labels <- c("1%","10%","100%")
    
  }
  
  #-----------------------------
  # INITIAL PLOT
  #-----------------------------
  
  p <- ggplot()
  
  if(nrow(other_data)>0){
    
    if(log_relative_prop){
      
      p <- p + geom_point(
        data=other_data,
        aes(x=x_num,
            y=y_num,
            size=log10(emission_rate_numeric+1),
            color=log_color),
        alpha=0.8)
      
    } else {
      
      p <- p + geom_point(
        data=other_data,
        aes(x=x_num,
            y=y_num,
            size=log10(emission_rate_numeric+1),
            color=relative_proportion),
        alpha=0.8)
      
    }
    
    p <- p +
      scale_size_continuous(
        name="Log Monoterpenes EF\n(µg·g⁻¹·h⁻¹)",
        range=c(2,10),
        breaks=reference_sizes,
        labels=bubble_labels
      )
    
  }
  
  if(nrow(isoprene_data)>0){
    
    p <- p + ggnewscale::new_scale("size")
    
    if(log_relative_prop){
      
      p <- p + geom_point(
        data=isoprene_data,
        aes(x=x_num,
            y=y_num,
            size=log10(emission_rate_numeric+1),
            color=log_color),
        alpha=0.8)
      
    } else {
      
      p <- p + geom_point(
        data=isoprene_data,
        aes(x=x_num,
            y=y_num,
            size=log10(emission_rate_numeric+1),
            color=relative_proportion),
        alpha=0.8)
      
    }
    
    p <- p +
      scale_size_continuous(
        name="Log Isoprene EF\n(µg·g⁻¹·h⁻¹)",
        range=c(2,8),
        breaks=reference_sizes_isoprene,
        labels=bubble_labels_isoprene
      )
    
  }
  
  #-----------------------------
  # COLOR SCALE
  #-----------------------------
  
  if(log_relative_prop){
    
    p <- p +
      scale_color_viridis_c(
        name="Log relative abundance (%)",
        option="inferno",
        direction=-1,
        breaks=color_breaks,
        labels=color_labels)
    
  } else {
    
    p <- p +
      scale_color_viridis_c(
        name="Relative proportion (%)",
        option="inferno",
        direction=-1,
        limits=c(0,100))
    
  }
  
  #-----------------------------
  # AXES
  #-----------------------------
  
  p <- p +
    scale_x_continuous(
      breaks=sort(unique(data_filtered$x_num)),
      labels=levels(data_filtered$species),
      expand=c(0,0)
    ) +
    scale_y_continuous(
      breaks=sort(unique(data_filtered$y_num)),
      labels=levels(data_filtered$compound),
      expand=c(0,0)
    )
  
  #-----------------------------
  # THEME
  #-----------------------------
  
  p <- p +
    theme_minimal() +
    theme(
      axis.text.x = element_text(angle=45,hjust=1,size=8),
      axis.text.y = element_text(size=7),
      legend.position = legend_position,
      plot.margin = margin(l=70,r=5,t=5,b=5),
      plot.title = element_text(hjust=0.5,face="bold",size=14)
    ) +
    labs(
      x="Species",
      y="Isoprenoid compounds",
      title=title
    )
  
  #-----------------------------
  # SAVE PLOT
  #-----------------------------
  
  if(save_plot){
    
    ggsave(
      filename=file_name,
      plot=p,
      device="png",
      width=width,
      height=height,
      units="mm",
      dpi=dpi,
      bg="white"
    )
    
  }
  
  return(p)
}

my_groups <- list(
  c("α-Terpinene", "α-Terpineol", "γ-Terpinene", "Limonene", "Linalool", "β-Phellandrene"),
  c("2,3,6-Trimethyl-1,5-heptadiene (ACI)", "Cosmene", "Isoterpinolene", 
    "Neoalloocimene", "Thuja-2,4(10)-diene", "alpha-Phellandrene", 
    "p-Cymenene", "p-Mentha-1,5,8-triene", "unknown monoterpene 1", 
    "3-Carene", "β-Pinene", "α-Pinene", "Camphene", "β-Ocimene", 
    "Isoprene", "Terpinolene", "p-Cymene", "Sabinene"),
  c("1,8-Cineole", "Myrcene", "Camphor")
)

group_labels <- c("10-11", "10-11", "10-12")




proportion_relative <- function(data) {
  # Transformer en format long
  data_long <- data %>%
    tidyr::pivot_longer(-compound, names_to = "species", values_to = "valeur") %>%
    dplyr::mutate(valeur = ifelse(is.na(valeur), 0, valeur))
  
  # Calculer somme par espèce
  data_long <- data_long %>%
    dplyr::group_by(species) %>%
    dplyr::mutate(
      total_espece = sum(valeur),
      proportion = ifelse(total_espece > 0, (valeur / total_espece) * 100, 0)
    ) %>%
    ungroup()
  
  # Remettre en format large
  data_wide <- data_long %>%
    dplyr::select(compound, species, proportion) %>%
    tidyr::pivot_wider(names_from = species, values_from = proportion, values_fill = 0)
  
  return(data_wide)
}

calculer_somme_terpenes <- function(df) {
  # Identifier les indices des terpènes
  indices <- which(df$compound %in% c("α-Terpinene", "alpha-Phellandrene", "Myrcene", "α-Terpineol", 
                                      "γ-Terpinene", "Limonene", "Linalool", 
                                      "β-Phellandrene", "β-Ocimene", "Sabinene"))
  
  # Calculer les sommes
  sommes_especes <- colSums(df[indices, -1], na.rm = TRUE)
  total <- sum(sommes_especes)
  
  # Retourner les résultats
  resultats <- data.frame(
    espece = c(names(sommes_especes), "TOTAL_GENERAL"),
    somme = c(as.numeric(sommes_especes), total)
  )
  
  return(resultats)
}


#  test <- create_bubble_heatmap_2 (table_heat_map_complete,
#                                                                         log_relative_prop = TRUE,
#                                                                         title = "",
#                                                                        bubble_sizes = c(1,10,50),
#                                                                         bubble_labels = NULL,
#                                                                         compound_order = my_compound_order,
#                                                                         species_order = my_species_order,
#                                                                        isoprene_compound = "Isoprene",
#                                                                         bubble_sizes_isoprene = c(1,10,500),
#                                                                         bubble_labels_isoprene = NULL,
#                                                                        legend_position = "bottom")
# print(test)

