identify_windows <- function(df,
                             metric_col,
                             threshold_1, CI_inf_1, CI_upp_1,
                             threshold_2 = NULL, CI_inf_2 = NULL, CI_upp_2 = NULL,
                             species_col = "sp_list") {
  
  # Vérifications
  if (!is.data.frame(df)) stop("df doit être un data.frame.")
  if (!metric_col %in% names(df)) stop("La colonne metric_col n'existe pas.")
  if (!species_col %in% names(df)) stop("La colonne species_col n'existe pas.")
  if (!is.numeric(df[[metric_col]])) stop("La colonne métrique doit être numérique.")
  
  # Extraire et trier
  metric <- df[[metric_col]]
  ordre <- order(metric)
  df_tri <- df[ordre, , drop = FALSE]
  metric_tri <- metric[ordre]
  n <- nrow(df_tri)
  if (n == 0) stop("Aucune fenêtre fournie.")
  
  # Fonctions auxiliaires pour les seuils
  find_before <- function(seuil) {
    indices <- which(metric_tri < seuil)
    if (length(indices) == 0) return(NA_integer_)
    max(indices)  # dernière fenêtre avant le seuil
  }
  find_after <- function(seuil) {
    indices <- which(metric_tri > seuil)
    if (length(indices) == 0) return(NA_integer_)
    min(indices)  # première fenêtre après le seuil
  }
  
  # Indices de base
  idx_min <- 1L
  idx_max <- n
  
  # Premier seuil
  idx_before_inf1 <- find_before(CI_inf_1)
  idx_after_upp1 <- find_after(CI_upp_1)
  
  # Fenêtre médiane entre min et before_inf1
  if (!is.na(idx_min) && !is.na(idx_before_inf1)) {
    idx_mid_min_inf1 <- floor((idx_min + idx_before_inf1) / 2)
  } else {
    idx_mid_min_inf1 <- NA_integer_
  }
  
  # Fenêtre médiane entre after_upp1 et max
  if (!is.na(idx_after_upp1) && !is.na(idx_max)) {
    idx_mid_upp1_max <- floor((idx_after_upp1 + idx_max) / 2)
  } else {
    idx_mid_upp1_max <- NA_integer_
  }
  
  indices <- c(idx_min, idx_mid_min_inf1, idx_before_inf1,
               idx_after_upp1, idx_mid_upp1_max, idx_max)
  
  # Second seuil éventuel
  if (!is.null(threshold_2)) {
    if (is.null(CI_inf_2) || is.null(CI_upp_2)) {
      stop("Si threshold_2 est fourni, CI_inf_2 et CI_upp_2 doivent l'être aussi.")
    }
    idx_before_inf2 <- find_before(CI_inf_2)
    idx_after_upp2 <- find_after(CI_upp_2)
    
    # Médiane entre after_upp1 et before_inf2
    if (!is.na(idx_after_upp1) && !is.na(idx_before_inf2)) {
      idx_mid_upp1_inf2 <- floor((idx_after_upp1 + idx_before_inf2) / 2)
    } else {
      idx_mid_upp1_inf2 <- NA_integer_
    }
    
    # Médiane entre after_upp2 et max
    if (!is.na(idx_after_upp2) && !is.na(idx_max)) {
      idx_mid_upp2_max <- floor((idx_after_upp2 + idx_max) / 2)
    } else {
      idx_mid_upp2_max <- NA_integer_
    }
    
    indices <- c(idx_min, idx_mid_min_inf1, idx_before_inf1,
                 idx_after_upp1, idx_mid_upp1_inf2,
                 idx_before_inf2, idx_after_upp2,
                 idx_mid_upp2_max, idx_max)
  }
  
  # Nettoyage des indices
  indices <- indices[!is.na(indices)]
  indices <- unique(indices)   # une même ligne n'apparaît qu'une fois
  indices <- sort(indices)
  
  if (length(indices) == 0) stop("Aucune fenêtre valide identifiée.")
  
  # Résultat
  result <- df_tri[indices, , drop = FALSE]
  
  
  return(result)
}

plot_convex_hulls_thresholds <- function(resultats, merged_data_3T, 
                                         sp_id_col = "rowname",
                                         axes = list(c("Comp.1", "Comp.2"),
                                                     c("Comp.1", "Comp.3"),
                                                     c("Comp.1", "Comp.4"),
                                                     c("Comp.2", "Comp.3"),
                                                     c("Comp.2", "Comp.4"),
                                                     c("Comp.3", "Comp.4")),
                                         colors = c("blue", "red"),
                                         alpha_poly = 0.3,
                                         show_window_points = FALSE,
                                         background_points = TRUE,
                                         windows_to_plot = NULL,
                                         window_colors = NULL,
                                         ncol = 2) {
  
  # Vérifications et chargement des packages
  if (!require(ggplot2)) stop("ggplot2 est requis")
  if (!require(dplyr)) stop("dplyr est requis")
  if (!require(patchwork)) stop("patchwork est requis (installer avec install.packages('patchwork'))")
  
  # Vérifier la présence de sp_list dans resultats
  if (!"sp_list" %in% names(resultats)) 
    stop("resultats doit contenir une colonne 'sp_list'")
  
  # Ajouter un identifiant de fenêtre basé sur l'ordre (tri croissant)
  resultats <- resultats %>% mutate(window_id = row_number())
  
  # Filtrer les fenêtres à représenter
  if (!is.null(windows_to_plot)) {
    if (!is.numeric(windows_to_plot)) stop("windows_to_plot doit être un vecteur numérique")
    resultats <- resultats %>% filter(window_id %in% windows_to_plot)
    if (nrow(resultats) == 0) stop("Aucune fenêtre sélectionnée après filtrage")
  }
  
  n_windows <- nrow(resultats)
  
  # Gestion des couleurs
  if (!is.null(window_colors)) {
    if (length(window_colors) != n_windows) {
      stop("La longueur de window_colors doit être égale au nombre de fenêtres sélectionnées (", n_windows, ")")
    }
    resultats <- resultats %>% mutate(group = as.character(window_id),
                                      group_color = window_colors)
    group_colors <- setNames(window_colors, as.character(resultats$window_id))
    use_individual_colors <- TRUE
  } else {
    group_size <- 3
    n_groups <- ceiling(n_windows / group_size)
    resultats <- resultats %>% 
      mutate(group = as.character(rep(1:n_groups, each = group_size, length.out = n_windows)))
    if (length(colors) < n_groups) {
      colors <- rep(colors, length.out = n_groups)
    }
    group_colors <- setNames(colors[1:n_groups], as.character(1:n_groups))
    use_individual_colors <- FALSE
  }
  
  # Vérifier que merged_data_3T contient les colonnes nécessaires
  required_cols <- unique(unlist(axes))
  if (!all(required_cols %in% names(merged_data_3T))) {
    stop("merged_data_3T ne contient pas toutes les colonnes d'axes requises")
  }
  if (!sp_id_col %in% names(merged_data_3T)) 
    stop(paste("La colonne", sp_id_col, "est absente de merged_data_3T"))
  
  # --- Préparer les données d'arrière‑plan (une seule fois, toutes paires) ---
  bg_points_by_pair <- list()
  if (background_points) {
    for (j in seq_along(axes)) {
      ax <- axes[[j]]
      xcol <- ax[1]; ycol <- ax[2]
      ok <- complete.cases(merged_data_3T[[xcol]], merged_data_3T[[ycol]])
      bg_pts <- data.frame(
        x = merged_data_3T[[xcol]][ok],
        y = merged_data_3T[[ycol]][ok]
      )
      bg_points_by_pair[[j]] <- bg_pts
    }
  }
  
  # --- Calcul des convex hulls et préparation des données par paire d'axes ---
  plots_list <- list()
  
  for (j in seq_along(axes)) {
    ax <- axes[[j]]
    xcol <- ax[1]
    ycol <- ax[2]
    pair_name <- paste(xcol, "vs", ycol)
    
    # Données d'arrière‑plan pour cette paire
    bg_pts <- if (background_points) bg_points_by_pair[[j]] else NULL
    
    # Polygones pour cette paire
    hull_list_pair <- list()
    points_list_pair <- list()
    
    for (i in 1:nrow(resultats)) {
      sp_list_i <- resultats$sp_list[[i]]
      if (is.list(sp_list_i)) sp_list_i <- unlist(sp_list_i)
      if (is.character(sp_list_i) && length(sp_list_i) == 1 && grepl(",", sp_list_i)) {
        sp_list_i <- trimws(unlist(strsplit(sp_list_i, ",")))
      }
      sp_data <- merged_data_3T[merged_data_3T[[sp_id_col]] %in% sp_list_i, , drop = FALSE]
      
      # Convex hull
      if (nrow(sp_data) >= 3) {
        coords <- sp_data[, c(xcol, ycol), drop = FALSE]
        coords <- na.omit(coords)
        if (nrow(coords) >= 3) {
          hull_idx <- chull(coords)
          hull_coords <- coords[hull_idx, ]
          hull_coords <- rbind(hull_coords, hull_coords[1, ])
          hull_df <- data.frame(
            group = resultats$group[i],   # déjà un caractère
            x = hull_coords[[xcol]],
            y = hull_coords[[ycol]]
          )
          hull_list_pair[[length(hull_list_pair) + 1]] <- hull_df
        }
      }
      
      # Points (optionnel)
      if (show_window_points && nrow(sp_data) > 0) {
        ok <- complete.cases(sp_data[[xcol]], sp_data[[ycol]])
        if (any(ok)) {
          pts <- data.frame(
            group = resultats$group[i],   # déjà un caractère
            x = sp_data[[xcol]][ok],
            y = sp_data[[ycol]][ok]
          )
          points_list_pair[[length(points_list_pair) + 1]] <- pts
        }
      }
    }
    
    if (length(hull_list_pair) == 0) {
      warning(paste("Aucun polygone pour la paire", pair_name, "- graphique ignoré"))
      next
    }
    
    hull_data_pair <- do.call(rbind, hull_list_pair)
    points_data_pair <- if (length(points_list_pair) > 0) do.call(rbind, points_list_pair) else NULL
    
    # Construction du graphique pour cette paire
    p_pair <- ggplot()
    
    if (background_points && !is.null(bg_pts)) {
      p_pair <- p_pair + geom_point(data = bg_pts, aes(x = x, y = y),
                                    color = "gray80", size = 0.5, alpha = 0.6)
    }
    
    p_pair <- p_pair + geom_polygon(data = hull_data_pair,
                                    aes(x = x, y = y, fill = group),
                                    alpha = alpha_poly, color = "black", linewidth = 0.3)
    
    if (show_window_points && !is.null(points_data_pair)) {
      p_pair <- p_pair + geom_point(data = points_data_pair,
                                    aes(x = x, y = y, color = group),
                                    size = 1.2, alpha = 0.8) +
        scale_color_manual(values = group_colors, guide = "none")
    }
    
    p_pair <- p_pair + 
      scale_fill_manual(values = group_colors, 
                        name = ifelse(use_individual_colors, "Fenêtre", "Groupe de fenêtres")) +
      labs(x = xcol, y = ycol) +
      theme_minimal() +
      theme(legend.position = "bottom")
    
    plots_list[[j]] <- p_pair
  }
  
  if (length(plots_list) == 0) 
    stop("Aucun graphique n'a pu être généré (vérifiez vos données)")
  
  # Assemblage avec patchwork
  combined_plot <- wrap_plots(plots_list, ncol = ncol) +
    plot_layout(guides = "collect") &
    theme(legend.position = "bottom")
  
  print(combined_plot)
  invisible(combined_plot)
}
plot_regularity <- function(resultats, merged_data_3T,
                            sp_id_col = "rowname",
                            axes = list(c("Comp.1", "Comp.2"),
                                        c("Comp.1", "Comp.3"),
                                        c("Comp.1", "Comp.4"),
                                        c("Comp.2", "Comp.3"),
                                        c("Comp.2", "Comp.4"),
                                        c("Comp.3", "Comp.4")),
                            windows_to_plot = NULL,
                            window_colors = NULL,
                            default_colors = c("blue", "red", "green", "orange", "purple"),
                            alpha_segment = 0.5,
                            alpha_point = 0.8,
                            size_point = 2,
                            size_segment = 0.5,
                            background_points = TRUE) {
  
  # Chargement des packages nécessaires
  if (!require(ggplot2)) stop("ggplot2 est requis")
  if (!require(dplyr)) stop("dplyr est requis")
  if (!require(tidyr)) stop("tidyr est requis pour les combinaisons de paires")
  
  # Vérifications initiales
  if (!"sp_list" %in% names(resultats)) 
    stop("resultats doit contenir une colonne 'sp_list'")
  
  # Ajouter un identifiant de fenêtre basé sur l'ordre (métrique croissante)
  resultats <- resultats %>% mutate(window_id = row_number())
  
  # Filtrer les fenêtres à représenter
  if (!is.null(windows_to_plot)) {
    if (!is.numeric(windows_to_plot)) stop("windows_to_plot doit être un vecteur numérique")
    resultats <- resultats %>% filter(window_id %in% windows_to_plot)
    if (nrow(resultats) == 0) stop("Aucune fenêtre sélectionnée après filtrage")
  }
  
  n_windows <- nrow(resultats)
  
  # Gestion des couleurs des fenêtres
  if (!is.null(window_colors)) {
    if (length(window_colors) != n_windows) {
      stop("La longueur de window_colors doit être égale au nombre de fenêtres sélectionnées (", n_windows, ")")
    }
    col_map <- setNames(window_colors, resultats$window_id)
  } else {
    # Palette par défaut étendue si nécessaire
    if (n_windows > length(default_colors)) {
      default_colors <- rep(default_colors, length.out = n_windows)
    }
    col_map <- setNames(default_colors[1:n_windows], resultats$window_id)
  }
  
  # Vérifier merged_data_3T
  required_cols <- unique(unlist(axes))
  if (!all(required_cols %in% names(merged_data_3T))) {
    stop("merged_data_3T ne contient pas toutes les colonnes d'axes requises")
  }
  if (!sp_id_col %in% names(merged_data_3T)) 
    stop(paste("La colonne", sp_id_col, "est absente de merged_data_3T"))
  
  # --- Préparer les données d'arrière‑plan (toutes espèces en gris) ---
  bg_points <- NULL
  if (background_points) {
    bg_list <- list()
    for (j in seq_along(axes)) {
      ax <- axes[[j]]
      xcol <- ax[1]; ycol <- ax[2]
      ok <- complete.cases(merged_data_3T[[xcol]], merged_data_3T[[ycol]])
      bg_pts <- data.frame(
        axes_pair = paste(xcol, "vs", ycol, sep = " "),
        x = merged_data_3T[[xcol]][ok],
        y = merged_data_3T[[ycol]][ok]
      )
      bg_list[[j]] <- bg_pts
    }
    bg_points <- do.call(rbind, bg_list)
  }
  
  # --- Fonction pour générer tous les segments d'un ensemble de points ---
  make_segments <- function(coords, window_id, axes_pair, color) {
    if (nrow(coords) < 2) return(NULL)
    # Toutes les combinaisons de paires
    pairs <- combn(seq_len(nrow(coords)), 2, simplify = FALSE)
    seg_df <- do.call(rbind, lapply(pairs, function(p) {
      data.frame(
        window_id = window_id,
        axes_pair = axes_pair,
        x = coords[p[1], 1],
        y = coords[p[1], 2],
        xend = coords[p[2], 1],
        yend = coords[p[2], 2],
        color = color
      )
    }))
    return(seg_df)
  }
  
  # --- Stocker les segments et les points ---
  all_segments <- list()
  all_points <- list()
  
  for (i in 1:nrow(resultats)) {
    wid <- resultats$window_id[i]
    col_wid <- col_map[as.character(wid)]
    # Extraire la liste des espèces
    sp_list_i <- resultats$sp_list[[i]]
    if (is.list(sp_list_i)) sp_list_i <- unlist(sp_list_i)
    if (is.character(sp_list_i) && length(sp_list_i) == 1 && grepl(",", sp_list_i)) {
      sp_list_i <- trimws(unlist(strsplit(sp_list_i, ",")))
    }
    # Récupérer les coordonnées pour ces espèces
    sp_data <- merged_data_3T[merged_data_3T[[sp_id_col]] %in% sp_list_i, , drop = FALSE]
    if (nrow(sp_data) == 0) next
    
    for (j in seq_along(axes)) {
      ax <- axes[[j]]
      xcol <- ax[1]; ycol <- ax[2]
      # Ne garder que les espèces avec coordonnées complètes
      ok <- complete.cases(sp_data[[xcol]], sp_data[[ycol]])
      if (!any(ok)) next
      coords <- sp_data[ok, c(xcol, ycol), drop = FALSE]
      colnames(coords) <- c("x", "y")
      
      # Segments
      segs <- make_segments(coords, wid, paste(xcol, "vs", ycol, sep = " "), col_wid)
      if (!is.null(segs)) all_segments[[length(all_segments)+1]] <- segs
      
      # Points
      pts <- data.frame(
        window_id = wid,
        axes_pair = paste(xcol, "vs", ycol, sep = " "),
        x = coords$x,
        y = coords$y,
        color = col_wid
      )
      all_points[[length(all_points)+1]] <- pts
    }
  }
  
  if (length(all_segments) == 0 && length(all_points) == 0)
    stop("Aucune donnée à afficher (vérifiez les espèces et leurs coordonnées).")
  
  segments_data <- if (length(all_segments) > 0) do.call(rbind, all_segments) else NULL
  points_data <- if (length(all_points) > 0) do.call(rbind, all_points) else NULL
  
  # --- Construction du graphique ---
  p <- ggplot()
  
  # Arrière‑plan : toutes les espèces (gris clair)
  if (background_points && !is.null(bg_points)) {
    p <- p + geom_point(data = bg_points, aes(x = x, y = y),
                        color = "gray85", size = 1, alpha = 0.5)
  }
  
  # Segments (arêtes)
  if (!is.null(segments_data)) {
    p <- p + geom_segment(data = segments_data,
                          aes(x = x, y = y, xend = xend, yend = yend, color = color),
                          alpha = alpha_segment, linewidth = size_segment)
  }
  
  # Points des espèces des fenêtres (avec la même couleur)
  if (!is.null(points_data)) {
    p <- p + geom_point(data = points_data,
                        aes(x = x, y = y, color = color),
                        size = size_point, alpha = alpha_point)
  }
  
  # Échelle des couleurs : une couleur par fenêtre (basée sur les valeurs uniques dans color)
  unique_colors <- unique(c(segments_data$color, points_data$color))
  names(unique_colors) <- unique_colors
  p <- p + scale_color_manual(values = unique_colors, name = "Fenêtre")
  
  # Facettes par combinaison d'axes
  p <- p + facet_wrap(~ axes_pair, scales = "free", ncol = 2) +
    labs(x = "", y = "") +
    theme_minimal() +
    theme(legend.position = "bottom",
          strip.background = element_rect(fill = "lightgray"))
  
  print(p)
  invisible(p)
}

plot_fdis <- function(resultats, merged_data_3T,
                      sp_id_col = "rowname",
                      axes = list(c("Comp.1", "Comp.2"),
                                  c("Comp.1", "Comp.3"),
                                  c("Comp.1", "Comp.4"),
                                  c("Comp.2", "Comp.3"),
                                  c("Comp.2", "Comp.4"),
                                  c("Comp.3", "Comp.4")),
                      windows_to_plot = NULL,
                      window_colors = NULL,
                      default_colors = c("blue", "red", "green", "orange", "purple"),
                      alpha_segment = 0.6,
                      alpha_point = 0.8,
                      size_point = 2,
                      size_segment = 0.5,
                      size_centroid = 4,
                      background_points = TRUE) {
  
  # Chargement des packages
  if (!require(ggplot2)) stop("ggplot2 est requis")
  if (!require(dplyr)) stop("dplyr est requis")
  
  # Vérifications
  if (!"sp_list" %in% names(resultats)) 
    stop("resultats doit contenir une colonne 'sp_list'")
  
  # Identifiant de fenêtre basé sur l'ordre croissant de métrique
  resultats <- resultats %>% mutate(window_id = row_number())
  
  # Filtrer les fenêtres à représenter
  if (!is.null(windows_to_plot)) {
    if (!is.numeric(windows_to_plot)) stop("windows_to_plot doit être numérique")
    resultats <- resultats %>% filter(window_id %in% windows_to_plot)
    if (nrow(resultats) == 0) stop("Aucune fenêtre sélectionnée")
  }
  
  n_windows <- nrow(resultats)
  
  # Palette de couleurs
  if (!is.null(window_colors)) {
    if (length(window_colors) != n_windows)
      stop("window_colors doit avoir longueur ", n_windows)
    col_map <- setNames(window_colors, resultats$window_id)
  } else {
    if (n_windows > length(default_colors))
      default_colors <- rep(default_colors, length.out = n_windows)
    col_map <- setNames(default_colors[1:n_windows], resultats$window_id)
  }
  
  # Vérifier merged_data_3T
  required_cols <- unique(unlist(axes))
  if (!all(required_cols %in% names(merged_data_3T)))
    stop("merged_data_3T ne contient pas toutes les colonnes d'axes")
  if (!sp_id_col %in% names(merged_data_3T))
    stop(paste("Colonne", sp_id_col, "absente"))
  
  # --- Arrière‑plan (toutes espèces) ---
  bg_points <- NULL
  if (background_points) {
    bg_list <- list()
    for (j in seq_along(axes)) {
      ax <- axes[[j]]
      xcol <- ax[1]; ycol <- ax[2]
      ok <- complete.cases(merged_data_3T[[xcol]], merged_data_3T[[ycol]])
      bg_pts <- data.frame(
        axes_pair = paste(xcol, "vs", ycol, sep = " "),
        x = merged_data_3T[[xcol]][ok],
        y = merged_data_3T[[ycol]][ok]
      )
      bg_list[[j]] <- bg_pts
    }
    bg_points <- do.call(rbind, bg_list)
  }
  
  # --- Fonction pour calculer centroïde d'un ensemble de points ---
  compute_centroid <- function(coords) {
    if (nrow(coords) == 0) return(NULL)
    data.frame(x = mean(coords[,1]), y = mean(coords[,2]))
  }
  
  # --- Collecte des données : points, centroïdes, segments ---
  all_points <- list()
  all_centroids <- list()
  all_segments <- list()
  
  for (i in 1:nrow(resultats)) {
    wid <- resultats$window_id[i]
    col_wid <- col_map[as.character(wid)]
    
    # Liste des espèces
    sp_list_i <- resultats$sp_list[[i]]
    if (is.list(sp_list_i)) sp_list_i <- unlist(sp_list_i)
    if (is.character(sp_list_i) && length(sp_list_i) == 1 && grepl(",", sp_list_i)) {
      sp_list_i <- trimws(unlist(strsplit(sp_list_i, ",")))
    }
    sp_data <- merged_data_3T[merged_data_3T[[sp_id_col]] %in% sp_list_i, , drop = FALSE]
    if (nrow(sp_data) == 0) next
    
    for (j in seq_along(axes)) {
      ax <- axes[[j]]
      xcol <- ax[1]; ycol <- ax[2]
      # Garder coordonnées complètes
      ok <- complete.cases(sp_data[[xcol]], sp_data[[ycol]])
      if (!any(ok)) next
      coords <- sp_data[ok, c(xcol, ycol), drop = FALSE]
      colnames(coords) <- c("x", "y")
      if (nrow(coords) == 0) next
      
      # Points
      pts <- data.frame(
        window_id = wid,
        axes_pair = paste(xcol, "vs", ycol, sep = " "),
        x = coords$x,
        y = coords$y,
        color = col_wid
      )
      all_points[[length(all_points)+1]] <- pts
      
      # Centroïde
      cent <- compute_centroid(coords)
      if (!is.null(cent)) {
        cent_df <- data.frame(
          window_id = wid,
          axes_pair = paste(xcol, "vs", ycol, sep = " "),
          x = cent$x,
          y = cent$y,
          color = col_wid
        )
        all_centroids[[length(all_centroids)+1]] <- cent_df
        
        # Segments : chaque point → centroïde
        segs <- data.frame(
          window_id = wid,
          axes_pair = paste(xcol, "vs", ycol, sep = " "),
          x = coords$x,
          y = coords$y,
          xend = cent$x,
          yend = cent$y,
          color = col_wid
        )
        all_segments[[length(all_segments)+1]] <- segs
      }
    }
  }
  
  if (length(all_points) == 0)
    stop("Aucune donnée à afficher.")
  
  points_data <- do.call(rbind, all_points)
  centroids_data <- if (length(all_centroids) > 0) do.call(rbind, all_centroids) else NULL
  segments_data <- if (length(all_segments) > 0) do.call(rbind, all_segments) else NULL
  
  # --- Construction graphique ---
  p <- ggplot()
  
  # Arrière‑plan
  if (background_points && !is.null(bg_points)) {
    p <- p + geom_point(data = bg_points, aes(x = x, y = y),
                        color = "gray85", size = 1, alpha = 0.5)
  }
  
  # Segments (points → centroïde)
  if (!is.null(segments_data)) {
    p <- p + geom_segment(data = segments_data,
                          aes(x = x, y = y, xend = xend, yend = yend, color = color),
                          alpha = alpha_segment, linewidth = size_segment)
  }
  
  # Points des espèces
  p <- p + geom_point(data = points_data,
                      aes(x = x, y = y, color = color),
                      size = size_point, alpha = alpha_point)
  
  # Centroïdes (étoiles)
  if (!is.null(centroids_data)) {
    p <- p + geom_point(data = centroids_data,
                        aes(x = x, y = y, color = color),
                        shape = 8, size = size_centroid, stroke = 1.2)
  }
  
  # Échelle de couleur
  unique_colors <- unique(c(points_data$color, 
                            if (!is.null(centroids_data)) centroids_data$color else NULL,
                            if (!is.null(segments_data)) segments_data$color else NULL))
  names(unique_colors) <- unique_colors
  p <- p + scale_color_manual(values = unique_colors, name = "Fenêtre")
  
  # Facettes
  p <- p + facet_wrap(~ axes_pair, scales = "free", ncol = 2) +
    labs(x = "", y = "") +
    theme_minimal() +
    theme(legend.position = "bottom",
          strip.background = element_rect(fill = "lightgray"))
  
  print(p)
  invisible(p)
}


plot_originality <- function(resultats, merged_data_3T,
                             sp_id_col = "rowname",
                             axes = list(c("Comp.1", "Comp.2"),
                                         c("Comp.1", "Comp.3"),
                                         c("Comp.1", "Comp.4"),
                                         c("Comp.2", "Comp.3"),
                                         c("Comp.2", "Comp.4"),
                                         c("Comp.3", "Comp.4")),
                             windows_to_plot = NULL,
                             window_colors = NULL,
                             default_colors = c("blue", "red", "green", "orange", "purple"),
                             alpha_point = 0.8,
                             size_point = 2,
                             size_centroid_window = 4,
                             size_centroid_global = 5,
                             shape_centroid_window = 8,   # étoile
                             shape_centroid_global = 18,  # losange
                             color_centroid_global = "gold",
                             show_segments_to_centroid = TRUE,
                             alpha_segment = 0.5,
                             size_segment = 0.5,
                             background_points = TRUE) {
  
  # Chargement des packages
  if (!require(ggplot2)) stop("ggplot2 est requis")
  if (!require(dplyr)) stop("dplyr est requis")
  
  # Vérifications
  if (!"sp_list" %in% names(resultats)) 
    stop("resultats doit contenir une colonne 'sp_list'")
  
  # Identifiant de fenêtre basé sur l'ordre croissant de métrique
  resultats <- resultats %>% mutate(window_id = row_number())
  
  # Filtrer les fenêtres à représenter
  if (!is.null(windows_to_plot)) {
    if (!is.numeric(windows_to_plot)) stop("windows_to_plot doit être numérique")
    resultats <- resultats %>% filter(window_id %in% windows_to_plot)
    if (nrow(resultats) == 0) stop("Aucune fenêtre sélectionnée")
  }
  
  n_windows <- nrow(resultats)
  
  # Palette de couleurs pour les fenêtres
  if (!is.null(window_colors)) {
    if (length(window_colors) != n_windows)
      stop("window_colors doit avoir longueur ", n_windows)
    col_map <- setNames(window_colors, resultats$window_id)
  } else {
    if (n_windows > length(default_colors))
      default_colors <- rep(default_colors, length.out = n_windows)
    col_map <- setNames(default_colors[1:n_windows], resultats$window_id)
  }
  
  # Vérifier merged_data_3T
  required_cols <- unique(unlist(axes))
  if (!all(required_cols %in% names(merged_data_3T)))
    stop("merged_data_3T ne contient pas toutes les colonnes d'axes")
  if (!sp_id_col %in% names(merged_data_3T))
    stop(paste("Colonne", sp_id_col, "absente"))
  
  # --- Arrière‑plan : toutes espèces en gris clair ---
  bg_points <- NULL
  if (background_points) {
    bg_list <- list()
    for (j in seq_along(axes)) {
      ax <- axes[[j]]
      xcol <- ax[1]; ycol <- ax[2]
      ok <- complete.cases(merged_data_3T[[xcol]], merged_data_3T[[ycol]])
      bg_pts <- data.frame(
        axes_pair = paste(xcol, "vs", ycol, sep = " "),
        x = merged_data_3T[[xcol]][ok],
        y = merged_data_3T[[ycol]][ok]
      )
      bg_list[[j]] <- bg_pts
    }
    bg_points <- do.call(rbind, bg_list)
  }
  
  # --- Centroïde global (toutes espèces de merged_data_3T) ---
  # On calcule pour chaque combinaison d'axes le centroïde des valeurs non-NA
  global_centroids <- list()
  for (j in seq_along(axes)) {
    ax <- axes[[j]]
    xcol <- ax[1]; ycol <- ax[2]
    ok <- complete.cases(merged_data_3T[[xcol]], merged_data_3T[[ycol]])
    if (any(ok)) {
      cx <- mean(merged_data_3T[[xcol]][ok], na.rm = TRUE)
      cy <- mean(merged_data_3T[[ycol]][ok], na.rm = TRUE)
      global_centroids[[j]] <- data.frame(
        axes_pair = paste(xcol, "vs", ycol, sep = " "),
        x = cx,
        y = cy
      )
    }
  }
  global_centroids_data <- do.call(rbind, global_centroids)
  
  # --- Collecte des données par fenêtre : points, centroïdes, segments (optionnels) ---
  all_points <- list()
  all_centroids <- list()
  all_segments <- list()
  
  for (i in 1:nrow(resultats)) {
    wid <- resultats$window_id[i]
    col_wid <- col_map[as.character(wid)]
    
    # Récupérer la liste des espèces
    sp_list_i <- resultats$sp_list[[i]]
    if (is.list(sp_list_i)) sp_list_i <- unlist(sp_list_i)
    if (is.character(sp_list_i) && length(sp_list_i) == 1 && grepl(",", sp_list_i)) {
      sp_list_i <- trimws(unlist(strsplit(sp_list_i, ",")))
    }
    sp_data <- merged_data_3T[merged_data_3T[[sp_id_col]] %in% sp_list_i, , drop = FALSE]
    if (nrow(sp_data) == 0) next
    
    for (j in seq_along(axes)) {
      ax <- axes[[j]]
      xcol <- ax[1]; ycol <- ax[2]
      ok <- complete.cases(sp_data[[xcol]], sp_data[[ycol]])
      if (!any(ok)) next
      coords <- sp_data[ok, c(xcol, ycol), drop = FALSE]
      colnames(coords) <- c("x", "y")
      
      # Points
      pts <- data.frame(
        window_id = wid,
        axes_pair = paste(xcol, "vs", ycol, sep = " "),
        x = coords$x,
        y = coords$y,
        color = col_wid
      )
      all_points[[length(all_points)+1]] <- pts
      
      # Centroïde de la fenêtre
      cx <- mean(coords$x)
      cy <- mean(coords$y)
      cent_df <- data.frame(
        window_id = wid,
        axes_pair = paste(xcol, "vs", ycol, sep = " "),
        x = cx,
        y = cy,
        color = col_wid
      )
      all_centroids[[length(all_centroids)+1]] <- cent_df
      
      # Segments vers le centroïde de la fenêtre (si demandé)
      if (show_segments_to_centroid) {
        segs <- data.frame(
          window_id = wid,
          axes_pair = paste(xcol, "vs", ycol, sep = " "),
          x = coords$x,
          y = coords$y,
          xend = cx,
          yend = cy,
          color = col_wid
        )
        all_segments[[length(all_segments)+1]] <- segs
      }
    }
  }
  
  if (length(all_points) == 0)
    stop("Aucune donnée à afficher pour les fenêtres sélectionnées.")
  
  points_data <- do.call(rbind, all_points)
  centroids_data <- do.call(rbind, all_centroids)
  segments_data <- if (length(all_segments) > 0) do.call(rbind, all_segments) else NULL
  
  # --- Construction graphique ---
  p <- ggplot()
  
  # Arrière‑plan
  if (background_points && !is.null(bg_points)) {
    p <- p + geom_point(data = bg_points, aes(x = x, y = y),
                        color = "gray85", size = 1, alpha = 0.5)
  }
  
  # Segments (optionnels)
  if (!is.null(segments_data)) {
    p <- p + geom_segment(data = segments_data,
                          aes(x = x, y = y, xend = xend, yend = yend, color = color),
                          alpha = alpha_segment, linewidth = size_segment)
  }
  
  # Points des espèces des fenêtres
  p <- p + geom_point(data = points_data,
                      aes(x = x, y = y, color = color),
                      size = size_point, alpha = alpha_point)
  
  # Centroïdes des fenêtres (étoiles)
  p <- p + geom_point(data = centroids_data,
                      aes(x = x, y = y, color = color),
                      shape = shape_centroid_window, size = size_centroid_window, stroke = 1)
  
  # Centroïde global (jaune, forme différente)
  if (!is.null(global_centroids_data) && nrow(global_centroids_data) > 0) {
    p <- p + geom_point(data = global_centroids_data,
                        aes(x = x, y = y),
                        shape = shape_centroid_global, size = size_centroid_global,
                        color = color_centroid_global, fill = color_centroid_global)
  }
  
  # Échelle de couleur pour les fenêtres
  unique_colors <- unique(c(points_data$color, centroids_data$color))
  names(unique_colors) <- unique_colors
  p <- p + scale_color_manual(values = unique_colors, name = "Fenêtre")
  
  # Facettes
  p <- p + facet_wrap(~ axes_pair, scales = "free", ncol = 2) +
    labs(x = "", y = "") +
    theme_minimal() +
    theme(legend.position = "bottom",
          strip.background = element_rect(fill = "lightgray"))
  
  print(p)
  invisible(p)
}

