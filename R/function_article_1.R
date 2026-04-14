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
  invisible(data.frame(
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
  
  # 5. Boucle
  for (i in seq_along(start_indices)) {
    idx <- start_indices[i]:(start_indices[i] + size - 1)
    subset <- data_clean[idx, ]
    
    median_ord[i]  <- median(subset[[ordination]], na.rm = TRUE)
    median_var1[i] <- median(subset[[var1]], na.rm = TRUE)
    median_var2[i] <- median(subset[[var2]], na.rm = TRUE)
  }
  
  # 6. Graphiques
  # Graphique 1 : médiane isoprène
  plot(median_ord, median_var1, 
       pch = 16, col = "blue",
       xlab = paste("Médiane de", ordination),
       ylab = paste("Médiane de", var1),
       main = paste("Évolution de la médiane de", var1))
  abline(h = mean(median_var1), lty = 2, col = "red")
  
  # Graphique 2 : médiane monoterpènes
  plot(median_ord, median_var2, 
       pch = 16, col = "darkgreen",
       xlab = paste("Médiane de", ordination),
       ylab = paste("Médiane de", var2),
       main = paste("Évolution de la médiane de", var2))
  abline(h = mean(median_var2), lty = 2, col = "red")
  
  # 7. Retourne un data.frame (invisible)
  invisible(data.frame(
    median_ordination = median_ord,
    median_isoprene   = median_var1,
    median_monoterpenes = median_var2
  ))
}




moving_window_convex <- function(data, ordination = "Sum", 
                                        size = 25, step = 1, col_scores, col_traits){
  
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
  
  
  
  
  for (i in seq_along(start_indices)) {
    idx <- start_indices[i]:(start_indices[i] + size - 1)
    subset <- data_clean[idx, ]
    
    subset_coords <- as.matrix(
      subset[c(col_scores)]
    )
    #richesse = volume convex hull
    volumes[i] <- convhulln(subset_coords, options = "FA")$vol
    
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
    
    dimentionality[i] <-var(pca.traits$sdev^2)
    
    ##valeur ordination
    valeurs_ord <- subset[[ordination]] 
    ordinations[i] <- median(valeurs_ord)
    
  }
  
  results <- setNames(
    data.frame(volumes, var_dist, f_dis, originality, dimentionality, ordinations),
    c("Richness", "regularity", "functional_dispersion", "originality", "dimensionality", paste0("median_", ordination))
  )
   
  return(results) 
}


plot_moving_window_convex <- function(data) {
  
  if (!is.data.frame(data)) stop("L'objet 'data' doit être un data.frame")
  if (ncol(data) < 2) stop("Le data.frame doit contenir au moins 2 colonnes")
  
  x_var <- names(data)[ncol(data)]                     # dernière colonne = variable x
  y_vars <- names(data)[1:min(5, ncol(data)-1)]        # métriques (max 5)
  
  # Liste pour stocker les graphiques
  plots <- list()
  
  for (y_var in y_vars) {
    p <- ggplot(data, aes(x = .data[[x_var]], y = .data[[y_var]])) +
      geom_point(size = 2, alpha = 0.7, color = "steelblue") +
      labs(
        title = paste("Évolution de", y_var, "en fonction de", x_var),
        x = x_var,
        y = y_var
      ) +
      theme_minimal() +
      theme(
        panel.grid.minor = element_blank(),
        plot.title = element_text(hjust = 0.5)
      )
    
    plots[[y_var]] <- p   # stocker avec le nom de la métrique
  }
  
  # Retourner la liste (l'utilisateur pourra afficher ou sauvegarder chaque plot)
  return(plots)
}
