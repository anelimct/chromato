simulate_spearman <- function(data, trait1 = "SLA", trait2 = "HeightMax", 
                              sample_sizes = 10:35, n_iter = 1000) {
  # Garder uniquement les espèces (lignes) sans valeurs manquantes pour les deux traits
  data_clean <- na.omit(data[, c(trait1, trait2)])
  
  # Initialiser une liste pour stocker les résultats
  results <- list()
  
  for (n in sample_sizes) {
    rho_vals <- numeric(n_iter)  # vecteur pour les 1000 rho
    
    for (i in 1:n_iter) {
      # Échantillonner n espèces sans remise
      idx <- sample(1:nrow(data_clean), size = n, replace = FALSE)
      subset <- data_clean[idx, ]
      
      # Calculer la corrélation de Spearman
      cor_test <- cor.test(subset[, trait1], subset[, trait2], method = "spearman")
      rho_vals[i] <- cor_test$estimate
    }
    
    # Stocker les résultats pour cette taille d'échantillon
    results[[as.character(n)]] <- data.frame(
      sample_size = n,
      mean_rho = mean(rho_vals),
      sd_rho = sd(rho_vals),
      all_rho = I(list(rho_vals)),  # I() pour conserver la liste dans le data.frame
      stringsAsFactors = FALSE
    )
  }
  
  # Combiner toutes les lignes en un seul data.frame
  do.call(rbind, results)
}


plot_spearman_sensibility <- function(resultats, title = NULL) {
  # Vérification de la structure des données
  if (!all(c("sample_size", "all_rho") %in% colnames(resultats))) {
    stop("Le data.frame 'resultats' doit contenir les colonnes 'sample_size' et 'all_rho'.")
  }
  
  # Calcul des quantiles et de la moyenne pour chaque taille d'échantillon
  percentiles <- do.call(rbind, lapply(1:nrow(resultats), function(i) {
    n <- resultats$sample_size[i]
    rho_vals <- resultats$all_rho[[i]]
    data.frame(
      sample_size = n,
      q02.5 = quantile(rho_vals, 0.025, na.rm = TRUE),
      q25   = quantile(rho_vals, 0.25, na.rm = TRUE),
      q50   = quantile(rho_vals, 0.5, na.rm = TRUE),
      q75   = quantile(rho_vals, 0.75, na.rm = TRUE),
      q97.5 = quantile(rho_vals, 0.975, na.rm = TRUE),
      mean_rho = mean(rho_vals, na.rm = TRUE)
    )
  }))
  
  # Titre par défaut si non spécifié
  if (is.null(title)) {
    title <- "Sensibilité de la corrélation de Spearman à la taille d'échantillon"
  }
  
  # Construction du graphique
  p <- ggplot(percentiles, aes(x = sample_size)) +
    geom_ribbon(aes(ymin = q25, ymax = q75), fill = "grey70", alpha = 0.5) +
    geom_ribbon(aes(ymin = q02.5, ymax = q97.5), fill = "grey90", alpha = 0.5) +
    geom_line(aes(y = q50), color = "black", size = 0.8, linetype = "solid") +
    geom_line(aes(y = mean_rho), color = "red", size = 0.8, linetype = "dashed") +
    geom_point(aes(y = mean_rho), color = "red", size = 1.5) +
    labs(
      title = title,
      x = "Taille d'échantillon (nombre d'espèces)",
      y = expression("ρ de Spearman"),
      caption = "1000 tirages avec remise. Zones grises : percentiles (2.5%-97.5% et 25%-75%). Ligne rouge : moyenne."
    ) +
    theme_minimal()
  
  return(p)
}

moving_window_spearman <- function(data, trait1 = "SLA", trait2 = "HeightMax",
                                   size = 20, step = 1, resultat,
                                   ordination = "Sum", decreasing = FALSE,
                                   plot_scatter_windows = FALSE,
                                   max_windows = NULL) {
  
  # 1. Ordonner selon la variable d'ordination
  data <- data[order(data[[ordination]], decreasing = decreasing), ]
  
  # 2. Nettoyage des données
  data_clean <- na.omit(data[, c(trait1, trait2, ordination)])
  
  n_total <- nrow(data_clean)
  if (size > n_total) stop("La taille de la fenêtre dépasse le nombre d'espèces disponibles.")
  
  # 3. Fenêtres
  start_indices <- seq(1, n_total - size + 1, by = step)
  n_windows <- length(start_indices)
  
  # 4. Initialisation
  rho_obs <- numeric(n_windows)
  median_ord <- numeric(n_windows)
  
  mean_trait1 <- numeric(n_windows)
  min_trait1  <- numeric(n_windows)
  max_trait1  <- numeric(n_windows)
  sd_trait1   <- numeric(n_windows)
  
  mean_trait2 <- numeric(n_windows)
  min_trait2  <- numeric(n_windows)
  max_trait2  <- numeric(n_windows)
  sd_trait2   <- numeric(n_windows)
  
  mean_global_trait1 <- mean(data_clean[[trait1]], na.rm = TRUE)
  mean_global_trait2 <- mean(data_clean[[trait2]], na.rm = TRUE)
  
  if (plot_scatter_windows) windows_data <- list()
  
  # 5. Boucle
  for (i in seq_along(start_indices)) {
    idx <- start_indices[i]:(start_indices[i] + size - 1)
    subset <- data_clean[idx, ]
    
    rho_obs[i] <- cor(subset[[trait1]], subset[[trait2]], method = "spearman")
    median_ord[i] <- median(subset[[ordination]])
    
    # Trait 1
    mean_trait1[i] <- mean(subset[[trait1]], na.rm = TRUE)
    min_trait1[i]  <- min(subset[[trait1]], na.rm = TRUE)
    max_trait1[i]  <- max(subset[[trait1]], na.rm = TRUE)
    sd_trait1[i]   <- sd(subset[[trait1]], na.rm = TRUE)
    
    # Trait 2
    mean_trait2[i] <- mean(subset[[trait2]], na.rm = TRUE)
    min_trait2[i]  <- min(subset[[trait2]], na.rm = TRUE)
    max_trait2[i]  <- max(subset[[trait2]], na.rm = TRUE)
    sd_trait2[i]   <- sd(subset[[trait2]], na.rm = TRUE)
    
    if (plot_scatter_windows) {
      subset$window_id <- i
      windows_data[[i]] <- subset
    }
  }
  
  # 6. Références simulées
  res_row <- resultat[resultat$sample_size == size, ]
  
  if (nrow(res_row) == 0) {
    warning("Pas de référence simulée")
    mean_rho <- NA
    sd_rho <- NA
  } else {
    mean_rho <- res_row$mean_rho
    sd_rho   <- res_row$sd_rho
  }
  
  # 7. SES
  if (!is.na(mean_rho) && !is.na(sd_rho) && sd_rho > 0) {
    ses <- (rho_obs - mean_rho) / sd_rho
  } else {
    ses <- rep(NA, n_windows)
  }
  
  # -------- GRAPH 1 : Corrélation --------
  plot(median_ord, rho_obs, pch = 16, col = "blue",
       xlab = paste("Médiane de", ordination),
       ylab = expression("ρ de Spearman"),
       main = paste("Corrélation glissante :", trait1, "~", trait2))
  abline(h = 0, lty = 2, col = "gray")
  
  if (!is.na(mean_rho)) {
    abline(h = mean_rho, col = "red", lwd = 2)
    abline(h = mean_rho + 1.96 * sd_rho, col = "red", lty = 2)
    abline(h = mean_rho - 1.96 * sd_rho, col = "red", lty = 2)
  }
  
  # Capturer le premier graphique
  first_plot <- recordPlot()
  
  # -------- GRAPH 2 : SES --------
  plot(median_ord, ses, pch = 16, col = "blue",
       xlab = paste("Médiane de", ordination),
       ylab = "SES",
       main = "SES de la corrélation")
  abline(h = 0, lty = 2, col = "gray")
  
  # -------- GRAPH 3 : Scatter (optionnel) --------
  if (plot_scatter_windows) {
    if (!is.null(max_windows)) windows_data <- windows_data[1:max_windows]
    
    for (i in seq_along(windows_data)) {
      sub <- windows_data[[i]]
      plot(sub[[trait1]], sub[[trait2]],
           main = paste("Fenêtre", i),
           xlab = trait1, ylab = trait2,
           pch = 16, col = "darkgreen")
    }
  }
  
  # -------- GRAPH 4 : Mean + min/max (trait1) --------
  plot(median_ord, mean_trait1, pch = 16, col = "blue",
       ylim = range(c(min_trait1, max_trait1)),
       xlab = paste("Médiane de", ordination),
       ylab = trait1,
       main = paste("Mean + range -", trait1))
  
  arrows(median_ord, min_trait1, median_ord, max_trait1,
         code = 3, angle = 90, length = 0.05, col = "gray")
  abline(h = mean_global_trait1, col = "red", lty = 2)
  
  # -------- GRAPH 5 : Mean + min/max (trait2) --------
  plot(median_ord, mean_trait2, pch = 16, col = "blue",
       ylim = range(c(min_trait2, max_trait2)),
       xlab = paste("Médiane de", ordination),
       ylab = trait2,
       main = paste("Mean + range -", trait2))
  
  arrows(median_ord, min_trait2, median_ord, max_trait2,
         code = 3, angle = 90, length = 0.05, col = "gray")
  abline(h = mean_global_trait2, col = "red", lty = 2)
  
  # -------- GRAPH 6 : SD trait1 --------
  plot(median_ord, sd_trait1, pch = 16, col = "purple",
       xlab = paste("Médiane de", ordination),
       ylab = paste("SD de", trait1),
       main = paste("Dispersion -", trait1))
  abline(h = mean(sd_trait1), col = "red", lty = 2)
  
  # -------- GRAPH 7 : SD trait2 --------
  plot(median_ord, sd_trait2, pch = 16, col = "purple",
       xlab = paste("Médiane de", ordination),
       ylab = paste("SD de", trait2),
       main = paste("Dispersion -", trait2))
  abline(h = mean(sd_trait2), col = "red", lty = 2)
  
  # -------- OUTPUT --------
    invisible(list(
      data = data.frame(
        median_ord = median_ord,
        rho_obs = rho_obs,
        ses = ses,
        mean_trait1 = mean_trait1,
        min_trait1 = min_trait1,
        max_trait1 = max_trait1,
        sd_trait1 = sd_trait1,
        mean_trait2 = mean_trait2,
        min_trait2 = min_trait2,
        max_trait2 = max_trait2,
        sd_trait2 = sd_trait2
      ),
      first_plot = first_plot
    ))
}

moving_window_medians <- function(data, 
                                  ordination = "Sum", 
                                  var1 = "isoprene", 
                                  var2 = "monoterpenes",
                                  size = 20, 
                                  step = 1) {
  
  # 1. Ordonner selon la variable d'ordination
  data <- data[order(data[[ordination]]), ]
  
  # 2. Garder uniquement les lignes complètes pour les trois colonnes
  data_clean <- na.omit(data[, c(ordination, var1, var2)])
  
  n_total <- nrow(data_clean)
  if (size > n_total) stop("La taille de la fenêtre dépasse le nombre d'observations.")
  
  # 3. Indices de début des fenêtres
  start_indices <- seq(1, n_total - size + 1, by = step)
  n_windows <- length(start_indices)
  
  # 4. Initialisation
  median_ord <- numeric(n_windows)
  median_var1 <- numeric(n_windows)
  median_var2 <- numeric(n_windows)
  sd_var1 <- numeric(n_windows)
  sd_var2 <- numeric(n_windows)
  
  # 5. Boucle
  for (i in seq_along(start_indices)) {
    idx <- start_indices[i]:(start_indices[i] + size - 1)
    subset <- data_clean[idx, ]
    
    median_ord[i]  <- median(subset[[ordination]], na.rm = TRUE)
    median_var1[i] <- median(subset[[var1]], na.rm = TRUE)
    median_var2[i] <- median(subset[[var2]], na.rm = TRUE)
    sd_var1[i]     <- sd(subset[[var1]], na.rm = TRUE)
    sd_var2[i]     <- sd(subset[[var2]], na.rm = TRUE)
  }
  
  # 6. Graphiques avec barres d'erreur (médiane ± écart-type)
  
  # Graphique 1 : médiane isoprène ± sd
  ylim1 <- range(c(median_var1 - sd_var1, median_var1 + sd_var1), na.rm = TRUE)
  plot(median_ord, median_var1, 
       pch = 16, col = "blue",
       xlab = paste("Médiane de", ordination),
       ylab = paste("Médiane de", var1),
       main = paste("Évolution de la médiane de", var1, "± écart-type"),
       ylim = ylim1)
  # Barres d'erreur : médiane ± écart-type
  segments(x0 = median_ord, y0 = median_var1 - sd_var1, 
           x1 = median_ord, y1 = median_var1 + sd_var1, 
           col = "gray60", lwd = 1.5)
  points(median_ord, median_var1, pch = 16, col = "blue")
  abline(h = mean(median_var1), lty = 2, col = "red")
  
  # Graphique 2 : médiane monoterpènes ± sd
  ylim2 <- range(c(median_var2 - sd_var2, median_var2 + sd_var2), na.rm = TRUE)
  plot(median_ord, median_var2, 
       pch = 16, col = "darkgreen",
       xlab = paste("Médiane de", ordination),
       ylab = paste("Médiane de", var2),
       main = paste("Évolution de la médiane de", var2, "± écart-type"),
       ylim = ylim2)
  segments(x0 = median_ord, y0 = median_var2 - sd_var2, 
           x1 = median_ord, y1 = median_var2 + sd_var2, 
           col = "gray60", lwd = 1.5)
  points(median_ord, median_var2, pch = 16, col = "darkgreen")
  abline(h = mean(median_var2), lty = 2, col = "red")
  
  # 7. Retourne un data.frame (invisible) incluant les écart-types
  invisible(data.frame(
    median_ordination = median_ord,
    median_isoprene   = median_var1,
    sd_isoprene       = sd_var1,
    median_monoterpenes = median_var2,
    sd_monoterpenes     = sd_var2
  ))
}




moving_window_convex <- function(data, ordination = "Sum", 
                                        size = 25, step = 1, col_scores, col_traits, col_species = "rowname"){
  
  # 1. Ordonner selon la variable d'ordination
  data <- data[order(data[[ordination]]), ]
  
  # 2. Garder uniquement les lignes avec la variable renseignée
  data_clean <- data[!is.na(data[[ordination]]), ]
  

  n_total <- nrow(data_clean)
  
  # 3. Indices de début des fenêtres
  start_indices <- seq(1, n_total - size + 1, by = step)
  n_windows <- length(start_indices)
  
  # Centroïde global
  coords <- as.matrix(data_clean[, c(col_scores)])
  centroid <- colMeans(coords)
  
  # Initialiser vecteur résultat
  volumes <- numeric(n_windows)
  var_dist <- numeric(n_windows)
  f_dis <- numeric(n_windows)
  originality <- numeric(n_windows)
  dimentionality <- numeric(n_windows)
  
  ordinations <- numeric(n_windows)
  
  # Liste pour stocker les espèces de chaque fenêtre
  species_lists <- vector("list", n_windows)
  
  
  for (i in seq_along(start_indices)) {
    idx <- start_indices[i]:(start_indices[i] + size - 1)
    subset <- data_clean[idx, ]
    
    
    species_lists[[i]] <- as.character(subset[[col_species]])
    
    
    subset_coords <- as.matrix(
      subset[c(col_scores)]
    )
    #richesse = volume convex hull
    volumes[i] <- geometry::convhulln(subset_coords, options = "FA")$vol
    
    #Régularité
    mat_dist <- dist(subset_coords)
    
    var_dist[i] <- var(as.vector(mat_dist))
    
    #diversity = functional dipersion (laliberté 2010)
    f_dis[i] <- fundiversity::fd_fdis(subset_coords)$FDis
    
    #originalité 
    centoids_points <- colMeans(subset_coords)
    originality[i] <- distance <- sqrt(sum((centoids_points - centroid )^2))
    
    #dimensionalité
    
    subset_traits <- as.data.frame (
      subset[c(col_traits)]
    )
    
    subset_traits_N <- subset_traits |> normaliser_dataframe()
    pca.traits<- princomp(subset_traits_N)
    
    dimentionality[i] <- var(pca.traits$sdev^2)
    
    ##valeur ordination
    valeurs_ord <- subset[[ordination]] 
    ordinations[i] <- median(valeurs_ord)
    
  }
  
  results <- setNames(
    data.frame(volumes, var_dist, f_dis, originality, dimentionality, ordinations, species_list = I(species_lists)),
    c("Richness", "regularity_inv", "functional_dispersion", "originality", "dimensionality_inv", paste0("median_", ordination), "sp_list")
  )
  
  results <- results |> 
    dplyr::mutate(
      regularity = 1 - (
        (regularity_inv - min(regularity_inv)) /
          (max(regularity_inv) - min(regularity_inv))
      )
    ) |> dplyr::mutate(dimensionality =  1 - (
                           (dimensionality_inv - min(dimensionality_inv)) /
                             (max(dimensionality_inv) - min(dimensionality_inv))
                           ))|>  dplyr::select(-regularity_inv, -dimensionality_inv)
  
  col_to_move <- paste0("median_", ordination)
  
  results <- results |>
    dplyr::relocate(all_of(col_to_move), .after = dplyr::last_col()) |>
    dplyr::relocate(sp_list, .after = dplyr::last_col())
  
   
  return(results) 
}

plot_moving_window_convex <- function(data, ncol = 2) {
  
  if (!is.data.frame(data)) stop("L'objet 'data' doit être un data.frame")
  if (ncol(data) < 2) stop("Le data.frame doit contenir au moins 2 colonnes")
  
  x_var <- names(data)[ncol(data) -1]                     # dernière colonne = variable x
  y_vars <- names(data)[1:min(5, ncol(data)-2)]        # métriques (max 5)
  
  # Liste pour stocker les graphiques
  plots <- list()
  
  for (y_var in y_vars) {
    p <- ggplot(data, aes(x = .data[[x_var]], y = .data[[y_var]])) +
      geom_point(size = 2, alpha = 0.7, color = "steelblue") +
      labs(
        title = y_var,  # Titre simplifié : juste le nom de la métrique
        x = x_var,
        y = y_var
      ) +
      theme_minimal() +
      theme(
        panel.grid.minor = element_blank(),
        plot.title = element_text(
          hjust = 0.5,
          background = element_rect(fill = "gray90", color = "gray70", 
                                    linetype = "solid", linewidth = 0.5),
          margin = margin(t = 2, b = 2, l = 5, r = 5)
        )
      )
    
    plots[[y_var]] <- p
  }
  
  # Affichage en grille
  if (length(plots) > 0) {
    n_plots <- length(plots)
    nrow_ <- ceiling(n_plots / ncol)
    
    # Vérifier si gridExtra est disponible
    if (!requireNamespace("gridExtra", quietly = TRUE)) {
      warning("Le package 'gridExtra' n'est pas installé. Affichage individuel des plots.")
      for (p in plots) print(p)
    } else {
      # Option : ouvrir une fenêtre plus grande (fonctionne dans certains environnements)
      if (interactive() && .Platform$OS.type != "windows") {
        try(dev.new(width = 12, height = 8), silent = TRUE)
      }
      gridExtra::grid.arrange(grobs = plots, ncol = ncol, nrow = nrow_)
    }
  }
  
  invisible(plots)
}


moving_window_categories <- function(data,
                                     ordination = "Sum",
                                     var1 = "isoprene",
                                     var2 = "monoterpenes",
                                     size = 20,
                                     step = 1,
                                     threshold_iso_emit = 1,
                                     threshold_mono_emit = 0.1) {
  
  # Vérification de ggplot2
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("La fonction nécessite le package 'ggplot2'. Installez-le avec install.packages('ggplot2').")
  }
  library(ggplot2)
  
  # Seuils pour les niveaux d'émission (définis par l'utilisateur dans sa question)
  seuil_iso_low <- 10
  seuil_iso_medium <- 30
  seuil_mono_low <- 2
  seuil_mono_medium <- 5.1
  
  # 1. Ordonner selon la variable d'ordination
  data <- data[order(data[[ordination]]), ]
  
  # 2. Garder uniquement les lignes complètes pour les trois colonnes
  data_clean <- na.omit(data[, c(ordination, var1, var2)])
  
  n_total <- nrow(data_clean)
  if (size > n_total) stop("La taille de la fenêtre dépasse le nombre d'observations.")
  
  # 3. Indices de début des fenêtres
  start_indices <- seq(1, n_total - size + 1, by = step)
  n_windows <- length(start_indices)
  
  # 4. Fonction de classification d'une observation
  classer <- function(iso, mono) {
    emit_iso <- !is.na(iso) && iso > threshold_iso_emit
    emit_mono <- !is.na(mono) && mono > threshold_mono_emit
    
    # Cas NE
    if (!emit_iso && !emit_mono) return("NE")
    
    # Niveaux pour iso et mono
    niveau_iso <- if (emit_iso) {
      if (iso < seuil_iso_low) "low"
      else if (iso <= seuil_iso_medium) "medium"
      else "high"
    } else NA
    
    niveau_mono <- if (emit_mono) {
      if (mono < seuil_mono_low) "low"
      else if (mono <= seuil_mono_medium) "medium"
      else "high"
    } else NA
    
    # Cas iso seulement
    if (emit_iso && !emit_mono) return(paste0("iso_", niveau_iso))
    # Cas mono seulement
    if (!emit_iso && emit_mono) return(paste0("mono_", niveau_mono))
    # Cas both
    if (emit_iso && emit_mono) return(paste0("both_", niveau_iso, "-", niveau_mono))
  }
  
  # Liste de toutes les sous-catégories possibles (ordre pour le graphique)
  all_subtypes <- c("NE",
                    "iso_low", "iso_medium", "iso_high",
                    "mono_low", "mono_medium", "mono_high",
                    "both_low-low", "both_low-medium", "both_low-high",
                    "both_medium-low", "both_medium-medium", "both_medium-high",
                    "both_high-low", "both_high-medium", "both_high-high")
  
  # Couleurs fournies par l'utilisateur
  subtype_colors <- c(
    "NE" = "#fce72e",
    "both_high-high" = "#c34f70",
    "both_high-medium" = "#eba5b5",
    "both_high-low" = "#f4c2cd",
    "both_medium-high" = "#eba5b5",      # approximation, non spécifiée
    "both_medium-medium" = "#eba5b5",
    "both_medium-low" = "#eba5b5",
    "both_low-high" = "#f4c2cd",
    "both_low-medium" = "#f4c2cd",
    "both_low-low" = "#f4c2cd",
    "iso_high" = "#59ac24",
    "iso_medium" = "#adda96",
    "iso_low" = "#adda96",               # non spécifiée, on prend medium
    "mono_high" = "#2c52c4",
    "mono_medium" = "#91a2e9",
    "mono_low" = "#b7c2f3"
  )
  # Compléter les couleurs manquantes pour les both non listés
  both_missing <- setdiff(grep("^both_", all_subtypes, value = TRUE), names(subtype_colors))
  for (b in both_missing) {
    if (grepl("high-", b)) subtype_colors[b] <- "#c34f70"
    else if (grepl("medium-", b)) subtype_colors[b] <- "#eba5b5"
    else if (grepl("low-", b)) subtype_colors[b] <- "#f4c2cd"
  }
  
  # 5. Initialisation
  median_ord <- numeric(n_windows)
  # Matrice des proportions : lignes = sous-catégories, colonnes = fenêtres
  prop_matrix <- matrix(0, nrow = length(all_subtypes), ncol = n_windows,
                        dimnames = list(all_subtypes, NULL))
  
  # 6. Boucle sur les fenêtres
  for (i in seq_len(n_windows)) {
    idx <- start_indices[i]:(start_indices[i] + size - 1)
    subset <- data_clean[idx, ]
    
    # Classification de chaque observation
    classes <- mapply(classer, subset[[var1]], subset[[var2]], SIMPLIFY = TRUE)
    
    # Tableau des fréquences pour toutes les catégories possibles
    tab <- table(factor(classes, levels = all_subtypes))
    
    # Proportion dans cette fenêtre
    prop_matrix[, i] <- as.vector(tab) / sum(tab)
    
    # Médiane de l'ordination
    median_ord[i] <- median(subset[[ordination]], na.rm = TRUE)
  }
  
  # 7. Création du data.frame pour ggplot2
  df_plot <- data.frame(
    window = rep(seq_len(n_windows), each = length(all_subtypes)),
    median_ordination = rep(median_ord, each = length(all_subtypes)),
    subtype = rep(all_subtypes, times = n_windows),
    proportion = as.vector(prop_matrix)
  )
  
  # Supprimer les lignes avec proportion nulle (optionnel, mais allège)
  df_plot <- df_plot[df_plot$proportion > 0, ]
  
  # 8. Graphique en barres empilées
  p <- ggplot(df_plot, aes(x = median_ordination, y = proportion, fill = subtype)) +
    geom_col(width = diff(range(median_ord)) / (n_windows * 1.2), 
             position = "stack", color = "black", size = 0.2) +
    scale_fill_manual(values = subtype_colors, name = "Sous-catégorie") +
    labs(x = paste("Médiane de", ordination),
         y = "Proportion",
         title = paste("Types d'émetteurs par fenêtre glissante\n(taille =", size, ", pas =", step, ")"),
         fill = "Catégorie") +
    theme_minimal() +
    theme(legend.position = "bottom",
          legend.text = element_text(size = 8),
          legend.key.size = unit(0.4, "cm"))
  
  print(p)
  
  # 9. Retour invisible des données
  invisible(list(proportions = prop_matrix,
                 median_ordination = median_ord,
                 data_plot = df_plot,
                 all_subtypes = all_subtypes))
}




extraire_especes_threshold <- function(sum_df, var_ordination, 
                                        lower, upper) {
  # Vérifications
  if (!is.data.frame(sum_df)) stop("sum_df doit être un data.frame")
  if (!var_ordination %in% colnames(sum_df)) 
    stop(paste("La variable", var_ordination, "n'est pas dans sum_df"))
  if (!"sp_list" %in% colnames(sum_df)) 
    stop("sum_df doit contenir une colonne 'sp_list' (liste d'espèces par fenêtre)")
  
  # 1. Sous-ensembles de fenêtres selon la valeur de la variable d'ordination
  low_windows <- sum_df[sum_df[[var_ordination]] < lower, ]
  high_windows <- sum_df[sum_df[[var_ordination]] > upper, ]
  mid_windows  <- sum_df[sum_df[[var_ordination]] >= lower & 
                           sum_df[[var_ordination]] <= upper, ]
  
  # 2. Fonction pour extraire les espèces uniques d'un data.frame de fenêtres
  get_unique_species <- function(df) {
    if (nrow(df) == 0) return(character(0))
    unique(unlist(df$sp_list))
  }
  
  species_low  <- get_unique_species(low_windows)
  species_high <- get_unique_species(high_windows)
  species_mid  <- get_unique_species(mid_windows)
  
  # 3. Espèces exclusives à chaque zone (présentes dans une seule)
  low_exclusive  <- setdiff(species_low,  union(species_high, species_mid))
  high_exclusive <- setdiff(species_high, union(species_low,  species_mid))
  
  # Retourner une liste nommée
  return(list(
    pre_seuil  = low_exclusive,
    seuil      = species_mid ,
    post_seuil = high_exclusive
  ))
}

moving_window_categories_2<- function(data,
                                                               ordination = "Sum",
                                                               var1 = "isoprene",
                                                               var2 = "monoterpenes",
                                                               size = 20,
                                                               step = 1,
                                                               threshold_iso_emit = 1,
                                                               threshold_mono_emit = 0.1,
                                                               plot_type = c("area", "bar"),
                                                               category_order = c("NE", "mono", "iso", "both"),
                                                               show_legend = FALSE) {
  
  plot_type <- match.arg(plot_type)
  
  # Vérification des packages
  if (plot_type == "area" && !requireNamespace("areaplot", quietly = TRUE)) {
    stop("Le mode 'area' nécessite le package 'areaplot'.")
  }
  if (plot_type == "bar" && !requireNamespace("ggplot2", quietly = TRUE)) {
    stop("Le mode 'bar' nécessite le package 'ggplot2'.")
  }
  
  # Seuils d'émission
  seuil_iso_low <- 10
  seuil_iso_medium <- 30
  seuil_mono_low <- 2
  seuil_mono_medium <- 5.1
  
  # Tri et nettoyage
  data <- data[order(data[[ordination]]), ]
  data_clean <- na.omit(data[, c(ordination, var1, var2)])
  n_total <- nrow(data_clean)
  if (size > n_total) stop("La taille de la fenêtre dépasse le nombre d'observations.")
  
  # Fenêtres glissantes
  start_indices <- seq(1, n_total - size + 1, by = step)
  n_windows <- length(start_indices)
  
  # Classification
  classer <- function(iso, mono) {
    emit_iso <- !is.na(iso) && iso > threshold_iso_emit
    emit_mono <- !is.na(mono) && mono > threshold_mono_emit
    
    if (!emit_iso && !emit_mono) return("NE")
    
    niveau_iso <- if (emit_iso) {
      if (iso < seuil_iso_low) "low"
      else if (iso <= seuil_iso_medium) "medium"
      else "high"
    } else NA
    
    niveau_mono <- if (emit_mono) {
      if (mono < seuil_mono_low) "low"
      else if (mono <= seuil_mono_medium) "medium"
      else "high"
    } else NA
    
    if (emit_iso && !emit_mono) return(paste0("iso_", niveau_iso))
    if (!emit_iso && emit_mono) return(paste0("mono_", niveau_mono))
    return(paste0("both_", niveau_iso, "-", niveau_mono))
  }
  
  # Noms des sous-catégories
  niveaux <- c("low", "medium", "high")
  all_mono <- paste0("mono_", niveaux)
  all_iso  <- paste0("iso_", niveaux)
  all_both <- paste0("both_", rep(niveaux, each = 3), "-", rep(niveaux, times = 3))
  
  # Ordre des catégories
  ordered_subtypes <- c()
  if ("NE" %in% category_order) ordered_subtypes <- c(ordered_subtypes, "NE")
  if ("mono" %in% category_order) ordered_subtypes <- c(ordered_subtypes, all_mono)
  if ("iso" %in% category_order)  ordered_subtypes <- c(ordered_subtypes, all_iso)
  if ("both" %in% category_order) ordered_subtypes <- c(ordered_subtypes, all_both)
  
  # Palette de couleurs : tous les both ont la même couleur
  col_palette <- c(
    "NE" = "#fce72e",
    "iso_low" = "#adda96",
    "iso_medium" = "#adda96",
    "iso_high" = "#59ac24",
    "mono_low" = "#b7c2f3",
    "mono_medium" = "#91a2e9",
    "mono_high" = "#2c52c4"
  )
  both_color <- "#c34f70"   # couleur unique pour tous les both
  for (b in all_both) col_palette[b] <- both_color
  subtype_colors <- col_palette[ordered_subtypes]
  
  # Calcul des proportions par fenêtre
  median_ord <- numeric(n_windows)
  prop_matrix <- matrix(0, nrow = length(ordered_subtypes), ncol = n_windows,
                        dimnames = list(ordered_subtypes, NULL))
  
  for (i in seq_len(n_windows)) {
    idx <- start_indices[i]:(start_indices[i] + size - 1)
    subset <- data_clean[idx, ]
    classes <- mapply(classer, subset[[var1]], subset[[var2]], SIMPLIFY = TRUE)
    tab <- table(factor(classes, levels = ordered_subtypes))
    prop_matrix[, i] <- tab / sum(tab)
    median_ord[i] <- median(subset[[ordination]], na.rm = TRUE)
  }
  
  # --- Correction des doublons sur l'axe des x ---
  # Agrégation des fenêtres ayant la même médiane (moyenne des proportions)
  df_temp <- data.frame(median = median_ord, t(prop_matrix))
  df_agg <- aggregate(. ~ median, data = df_temp, FUN = mean)
  x_vals <- df_agg$median
  y_mat <- as.matrix(df_agg[, -1])
  colnames(y_mat) <- ordered_subtypes  # garantir les noms
  
  # Tri par médiane croissante
  ord_idx <- order(x_vals)
  x_vals <- x_vals[ord_idx]
  y_mat <- y_mat[ord_idx, , drop = FALSE]
  
  # Vérification finale (plus aucun doublon normalement)
  if (any(duplicated(x_vals))) stop("Il reste des doublons dans les abscisses.")
  
  # Tracé
  if (plot_type == "area") {
    old_par <- par(mar = c(7, 4, 4, 2) + 0.1)
    on.exit(par(old_par))
    
    areaplot::areaplot(x = x_vals, y = y_mat,
                       prop = FALSE, rev = FALSE,
                       col = subtype_colors, border = NA,
                       xlab = paste("Médiane de", ordination),
                       ylab = "Proportion",
                       main = paste("Types d'émetteurs\n(taille =", size, ", pas =", step, ")"))
    
    if (show_legend) {
      non_zero <- colSums(y_mat) > 0
      legend_labels <- ordered_subtypes[non_zero]
      legend_cols <- subtype_colors[non_zero]
      ncol_leg <- min(6, length(legend_labels))
      legend("bottom", inset = c(0, -0.25), legend = legend_labels,
             fill = legend_cols, ncol = ncol_leg, bty = "n", cex = 0.8, xpd = NA)
    }
  } else { # barres avec ggplot2
    library(ggplot2)
    n_windows_agg <- nrow(y_mat)
    df_plot <- data.frame(
      median_ordination = rep(x_vals, each = length(ordered_subtypes)),
      subtype = rep(ordered_subtypes, times = n_windows_agg),
      proportion = as.vector(t(y_mat))
    )
    df_plot <- df_plot[df_plot$proportion > 0, ]
    df_plot$subtype <- factor(df_plot$subtype, levels = ordered_subtypes)
    
    p <- ggplot(df_plot, aes(x = median_ordination, y = proportion, fill = subtype)) +
      geom_col(width = diff(range(x_vals)) / (n_windows_agg * 1.2),
               position = "stack", color = "black", size = 0.2) +
      scale_fill_manual(values = subtype_colors, name = "Catégorie") +
      labs(x = paste("Médiane de", ordination), y = "Proportion",
           title = paste("Types d'émetteurs\n(taille =", size, ", pas =", step, ")")) +
      theme_minimal() +
      theme(legend.position = if(show_legend) "bottom" else "none",
            legend.text = element_text(size = 8),
            legend.key.size = unit(0.4, "cm"))
    print(p)
  }
  
  invisible(list(proportions = prop_matrix,
                 median_ordination = median_ord,
                 aggregated_x = x_vals,
                 aggregated_prop = y_mat,
                 ordered_subtypes = ordered_subtypes))
}
