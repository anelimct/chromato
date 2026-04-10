
 library(readxl)
 MMIR_FRAC <- read_excel("data/MIR_FAC_bon.xlsx")
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
 # file: analyse_mir_fac.R
 
 # file: analyse_mir_fac.R
 analyse_mir_fac <- function(compound_data, mmir_frac_data, k = 3,
                             save_plots = TRUE,
                             plot_dir = "figures/screening",
                             label_size = 3.2) {
   
   library(dplyr)
   library(tidyr)
   library(ggplot2)
   library(ggrepel)
   
   # --------------------------------------------------------------------
   # 1. Data prep
   # --------------------------------------------------------------------
   compound_data_renamed <- compound_data %>%
     rename(Compound = compound)
   
   merged_data <- left_join(compound_data_renamed, mmir_frac_data, by = "Compound") %>%
     filter(!class %in% c("Sesquiterpenes", "Oxygenated-sesquiterpenes"))
   
   mean_columns <- names(merged_data)[grepl("^mean_", names(merged_data))]
   
   merged_data <- merged_data %>%
     mutate(across(all_of(mean_columns), ~ as.numeric(as.character(.)))) %>%
     mutate(
       MIR = as.numeric(as.character(MIR)),
       FAC = as.numeric(as.character(FAC))
     )
   
   # --------------------------------------------------------------------
   # 2. OFP / SOAFP
   # --------------------------------------------------------------------
   long_data <- merged_data %>%
     select(Compound, MIR, FAC, all_of(mean_columns)) %>%
     pivot_longer(cols = all_of(mean_columns),
                  names_to = "Species_col",
                  values_to = "E") %>%
     mutate(Species = gsub("^mean_", "", Species_col)) %>%
     select(-Species_col)
   
   contributions <- long_data %>%
     mutate(
       OFP_contrib = ifelse(!is.na(E) & !is.na(MIR), E * MIR, 0),
       SOAFP_contrib = ifelse(!is.na(E) & !is.na(FAC), E * (FAC / 100), 0)
     ) %>%
     group_by(Species) %>%
     summarise(
       OFP = sum(OFP_contrib, na.rm = TRUE),
       SOAFP = sum(SOAFP_contrib, na.rm = TRUE),
       .groups = "drop"
     )
   
   final_results <- contributions
   
   # --------------------------------------------------------------------
   # 3. Clustering
   # --------------------------------------------------------------------
   cluster_data <- final_results %>%
     select(OFP, SOAFP) %>%
     na.omit()
   
   species_names <- final_results$Species[
     !is.na(final_results$OFP) & !is.na(final_results$SOAFP)
   ]
   
   scaled_data <- scale(cluster_data)
   
   set.seed(123)
   km <- kmeans(scaled_data, centers = k, nstart = 25)
   
   final_results$Cluster <- NA
   final_results$Cluster[match(species_names, final_results$Species)] <- km$cluster
   final_results$Cluster <- factor(final_results$Cluster)
   
   # --------------------------------------------------------------------
   # 4. Labels propres (clean simple)
   # --------------------------------------------------------------------
   full_names <- c(
     bpub = "Betula pubescens",
     qcre = "Quercus crenata",
     vagn = "Vitex agnus-castus",
     bpen = "Betula pendula",
     jthu = "Juniperus thurifera",
     fcar = "Ficus carica",
     salb = "Salix alba",
     scin = "Salix cinerea",
     rala = "Rhamnus alaternus",
     rcat = "Rhamnus catharica",
     spur = "Salix purpurea",
     sele = "Salix eleagnos",
     spen = "Salix pentandra",
     faln = "Frangula alnus",
     ptre = "Populus tremula"
   )
   
   final_results <- final_results %>%
     mutate(label_display = ifelse(
       Species %in% names(full_names),
       full_names[Species],
       NA_character_
     ))
   
   # --------------------------------------------------------------------
   # 5. Colors
   # --------------------------------------------------------------------
   cluster_colors <- setNames(
     c("#fdbc5d", "#ad61b6", "#d04446")[1:k],
     as.character(1:k)
   )
   
   # --------------------------------------------------------------------
   # 6. Plot (FIXED + jitter effect + clean repel)
   # --------------------------------------------------------------------
   p_final <- ggplot(final_results,
                     aes(x = OFP, y = SOAFP,
                         label = label_display,
                         color = Cluster)) +
     
     geom_point(size = 3.5, alpha = 0.9) +
     
     geom_text_repel(
       aes(color = Cluster),
       
       size = label_size,          # ↓ slightly smaller (your request)
       family = "Arial",
       fontface = "italic",
       
       # --- KEY FIXES FOR OVERLAP ---
       box.padding = 0.7,
       point.padding = 0.5,
       
       # jitter-like behaviour
       nudge_x = 0,
       nudge_y = 0,
       
       # stronger separation
       force = 5,
       force_pull = 0.6,
       
       # connector lines (colored correctly)
       segment.color = NA,
       segment.size = 0.6,
       
       # avoid stacking
       max.overlaps = Inf,
       show.legend = FALSE
     ) +
     
     scale_color_manual(values = cluster_colors, name = "Cluster") +
     
     labs(
       x = expression("OFP (" * µ * "g·g"^{-1}*"·h"^{-1} * ")"),
       y = expression("SOAFP (" * µ * "g·g"^{-1}*"·h"^{-1} * ")")
     ) +
     
     theme_minimal(base_size = 12) +
     theme(
       text = element_text(family = "Arial"),
       legend.title = element_text(face = "bold")
     )
   
   print(p_final)
   
   # --------------------------------------------------------------------
   # 7. Save
   # --------------------------------------------------------------------
   if (save_plots) {
     dir.create(plot_dir, showWarnings = FALSE, recursive = TRUE)
     ggsave(file.path(plot_dir, "clustering_OFP_SOAFP.png"),
            plot = p_final, width = 8, height = 6, dpi = 300)
   }
   
   invisible(list(
     results = final_results,
     clustering = km,
     plot = p_final
   ))
 }
 