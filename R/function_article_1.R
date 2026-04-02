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
                                   ordination = "Sum", decreasing = FALSE) {
  
  # 1. Ordonner selon la variable d'ordination
  data <- data[order(data[[ordination]], decreasing = decreasing), ]
  
  # 2. Garder uniquement les lignes avec les deux traits renseignés
  data_clean <- na.omit(data[, c(trait1, trait2, ordination)])
  
  n_total <- nrow(data_clean)
  if (size > n_total) stop("La taille de la fenêtre dépasse le nombre d'espèces disponibles.")
  
  # 3. Indices de début des fenêtres
  start_indices <- seq(1, n_total - size + 1, by = step) # list des indexes de position de départ de la moving window, de 1 audernier indice où une fenêtre de taille size peut commencer sans dépasser la longueur totale des données.
  
  # 4. Calcul des corrélations et médianes
  rho_obs <- numeric(length(start_indices)) #initialisation des vecteurs
  median_ord <- numeric(length(start_indices))
  
  for (i in seq_along(start_indices)) {
    idx <- start_indices[i]:(start_indices[i] + size - 1) # position de départ et de fin de la moving window
    subset <- data_clean[idx, ]
    rho_obs[i] <- cor(subset[[trait1]], subset[[trait2]], method = "spearman")
    median_ord[i] <- median(subset[[ordination]])
  }
  
  # 5. Récupérer les valeurs simulées pour cette taille de fenêtre, a partir du bootstrap
  res_row <- resultat[resultat$sample_size == size, ]
  if (nrow(res_row) == 0) {
    warning("La taille de fenêtre n'est pas dans 'resultat' ; les lignes de référence seront absentes.")
    mean_rho <- NA
    sd_rho <- NA
  } else {
    mean_rho <- res_row$mean_rho
    sd_rho   <- res_row$sd_rho
  }
  
  # 6. Calcul du SES si possible, sinon NA
  if (!is.na(mean_rho) && !is.na(sd_rho) && sd_rho > 0) {
    ses <- (rho_obs - mean_rho) / sd_rho
  } else {
    ses <- rep(NA, length(rho_obs))
  }
  
  
  # ---- Graphique 1 : ρ observé + références simulées ----
  plot(median_ord, rho_obs, type = "p", pch = 16, col = "blue",
       xlab = paste("Médiane de", ordination, "dans la fenêtre"),
       ylab = expression("ρ de Spearman"),
       main = paste("Corrélation glissante :", trait1, "~", trait2,
                    "\nTaille =", size, ", pas =", step))
  abline(h = 0, lty = 2, col = "gray")
  if (!is.na(mean_rho)) {
    abline(h = mean_rho, col = "red", lwd = 2)
    abline(h = mean_rho + 1.96 * sd_rho, col = "red", lty = 2)
    abline(h = mean_rho - 1.96 * sd_rho, col = "red", lty = 2)
    legend("topright",  inset = c(-0.2, 0), xpd = TRUE,legend = c("ρ observé", "moyenne simulée", "moyenne ± 1,96σ"),
           col = c("blue", "red", "red"), lty = c(NA, 1, 2), lwd = c(1, 2, 1),
           pch = c(16, NA, NA), bty = "n")
  } else {
    legend("topright", legend = "ρ observé", col = "blue", pch = 16, bty = "n")
  }
  
  # ---- Graphique 2 : SES ----
  plot(median_ord, ses, type = "p", pch = 16, col = "blue",
       xlab = paste("Médiane de", ordination, "dans la fenêtre"),
       ylab = "SES (effet standardisé)",
       main = paste("SES de la corrélation glissante\nTaille =", size, ", pas =", step))
  abline(h = 0, lty = 2, col = "gray")

  
  # Retour silencieux des données
  invisible(data.frame(median_ord = median_ord, rho_obs = rho_obs, ses = ses))
}
