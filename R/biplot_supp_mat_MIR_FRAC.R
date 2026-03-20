
 library(readxl)
 MMIR_FRAC <- read_excel("data/MMIR_FRAC.xlsx")
 targets::tar_load(compound_mean_spagg)


#' MIR/FAC analysis: sum calculation, visualization, and clustering
#'
#' @param compound_data Dataframe with compounds (must contain columns
#'   `compound`, `class`, and `mean_XXX` columns for each species).
#' @param mmir_frac_data Dataframe with MIR and FAC factors (columns
#'   `Compound`, `MIR`, `FAC`).
#' @param k Number of clusters for clustering (default = 3). Adjust after
#'   inspecting the elbow plot.
#' @param save_plots Logical; if TRUE, saves plots to `plot_dir`.
#' @param plot_dir Destination folder for plots (created if needed).
#'
#' @return A list containing:
#'   - `results` : final dataframe with species sums
#'   - `correlation` : correlation coefficient between Sum_FAC_E and Sum_MIR_E
#'   - `plots` : list of ggplot objects (scatter, elbow, cluster comparison)
#'
#' @examples
#' \dontrun{
#'   results <- analyse_mir_fac(compound_mean_spagg_to_keep, MMIR_FRAC, k = 3)
#' }
analyse_mir_fac <- function(compound_data, mmir_frac_data, k = 3,
                            save_plots = FALSE, plot_dir = "figures/screening") {
  
  # Load required packages
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(ggrepel)
  library(patchwork)
  
  # --------------------------------------------------------------------
  # 1. Data preparation
  # --------------------------------------------------------------------
  # Rename compound column for join
  compound_data_renamed <- compound_data %>%
    rename(Compound = compound)
  
  # Merge with MIR/FAC factors and remove sesquiterpenes
  merged_data <- left_join(compound_data_renamed, mmir_frac_data, by = "Compound") %>%
    filter(!class %in% c("Sesquiterpenes", "Oxygenated-sesquiterpenes"))
  
  # Identify mean columns per species
  mean_columns <- names(merged_data)[grepl("^mean_", names(merged_data))]
  
  # --------------------------------------------------------------------
  # 2. Calculate MIR*E and FAC*E sums per species
  # --------------------------------------------------------------------
  mir_results <- data.frame(Species = character(), Sum_MIR_E = numeric())
  fac_results <- data.frame(Species = character(), Sum_FAC_E = numeric())
  
  for (esp_col in mean_columns) {
    species <- gsub("^mean_", "", esp_col)
    sum_mir <- 0
    sum_fac <- 0
    
    for (i in 1:nrow(merged_data)) {
      e_val <- as.numeric(merged_data[[esp_col]][i])
      mir_val <- merged_data$MIR[i]
      fac_val <- merged_data$FAC[i]
      
      if (!is.na(e_val) && !is.na(mir_val)) {
        sum_mir <- sum_mir + (e_val * mir_val)
      }
      if (!is.na(e_val) && !is.na(fac_val)) {
        sum_fac <- sum_fac + (e_val * fac_val)
      }
    }
    
    mir_results <- rbind(mir_results, data.frame(Species = species, Sum_MIR_E = sum_mir))
    fac_results <- rbind(fac_results, data.frame(Species = species, Sum_FAC_E = sum_fac))
  }
  
  final_results <- full_join(mir_results, fac_results, by = "Species")
  
  # --------------------------------------------------------------------
  # 3. Basic scatter plot and correlation
  # --------------------------------------------------------------------
  p_scatter <- ggplot(final_results, aes(x = Sum_FAC_E, y = Sum_MIR_E, label = Species)) +
    geom_point(size = 3, alpha = 0.7, color = "steelblue") +
    geom_text_repel(size = 3, max.overlaps = 20) +
    labs(
      title = "Sum(MIR·E) vs Sum(FAC·E)",
      x = "Sum(FAC·E)",
      y = "Sum(MIR·E)"
    ) +
    theme_minimal() +
    theme(
      plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
      axis.title = element_text(size = 12),
      axis.text = element_text(size = 10)
    )
  
  correlation <- cor(final_results$Sum_FAC_E, final_results$Sum_MIR_E,
                     use = "complete.obs")
  cat("\nCorrelation between Sum_FAC_E and Sum_MIR_E:", round(correlation, 3), "\n")
  
  # --------------------------------------------------------------------
  # 4. Clustering
  # --------------------------------------------------------------------
  # Prepare data (remove NAs and scale)
  cluster_data <- final_results %>%
    select(Sum_FAC_E, Sum_MIR_E) %>%
    na.omit()
  
  species_names <- final_results$Species[
    !is.na(final_results$Sum_FAC_E) & !is.na(final_results$Sum_MIR_E)
  ]
  
  scaled_data <- scale(cluster_data)
  
  # Elbow method (compute WSS for k = 1..10)
  wss <- sapply(1:10, function(kk) {
    kmeans(scaled_data, centers = kk, nstart = 25)$tot.withinss
  })
  
  df_elbow <- data.frame(k = 1:10, wss = wss)
  
  p_elbow <- ggplot(df_elbow, aes(x = k, y = wss)) +
    geom_line() +
    geom_point(size = 2) +
    geom_vline(xintercept = k, linetype = "dashed", color = "red") +
    labs(
      title = "Elbow method for k‑means",
      x = "Number of clusters k",
      y = "Within-cluster sum of squares (WSS)"
    ) +
    theme_minimal()
  
  # Hierarchical clustering (Ward)
  dist_matrix <- dist(scaled_data, method = "euclidean")
  hc <- hclust(dist_matrix, method = "ward.D")
  hc_clusters <- cutree(hc, k = k)
  
  # K‑means
  set.seed(123)
  km <- kmeans(scaled_data, centers = k, nstart = 25)
  km_clusters <- km$cluster
  
  # Add clusters to final dataframe
  final_results$Cluster_hc <- NA
  final_results$Cluster_hc[match(species_names, final_results$Species)] <- hc_clusters
  final_results$Cluster_hc <- factor(final_results$Cluster_hc)
  
  final_results$Cluster_km <- NA
  final_results$Cluster_km[match(species_names, final_results$Species)] <- km_clusters
  final_results$Cluster_km <- factor(final_results$Cluster_km)
  
  # Comparative cluster plots
  p_hc <- ggplot(final_results, aes(x = Sum_FAC_E, y = Sum_MIR_E,
                                    label = Species, color = Cluster_hc)) +
    geom_point(size = 3, alpha = 0.7) +
    geom_text_repel(size = 3, max.overlaps = 20, show.legend = FALSE) +
    scale_color_brewer(palette = "Set1") +
    labs(title = paste("Hierarchical (Ward) – k =", k),
         x = "Sum(FAC·E)", y = "Sum(MIR·E)", color = "Cluster") +
    theme_minimal()
  
  p_km <- ggplot(final_results, aes(x = Sum_FAC_E, y = Sum_MIR_E,
                                    label = Species, color = Cluster_km)) +
    geom_point(size = 3, alpha = 0.7) +
    geom_text_repel(size = 3, max.overlaps = 20, show.legend = FALSE) +
    scale_color_brewer(palette = "Set1") +
    labs(title = paste("K‑means – k =", k),
         x = "Sum(FAC·E)", y = "Sum(MIR·E)", color = "Cluster") +
    theme_minimal()
  
  p_cluster_compare <- p_hc + p_km +
    plot_annotation(title = "Comparison of clustering methods")
  
  # --------------------------------------------------------------------
  # 5. Optional plot saving
  # --------------------------------------------------------------------
  if (save_plots) {
    dir.create(plot_dir, showWarnings = FALSE, recursive = TRUE)
    
    ggsave(file.path(plot_dir, "scatter_mir_fac.png"),
           plot = p_scatter, width = 8, height = 6, dpi = 300)
    ggsave(file.path(plot_dir, "elbow_plot.png"),
           plot = p_elbow, width = 6, height = 4, dpi = 300)
    ggsave(file.path(plot_dir, "cluster_comparison.png"),
           plot = p_cluster_compare, width = 12, height = 6, dpi = 300)
    
    message("Plots saved to: ", plot_dir)
  }
  
  # --------------------------------------------------------------------
  # 6. Return results
  # --------------------------------------------------------------------
  invisible(list(
    results = final_results,
    correlation = correlation,
    plots = list(
      scatter = p_scatter,
      elbow = p_elbow,
      cluster_compare = p_cluster_compare
    )
  ))
}
