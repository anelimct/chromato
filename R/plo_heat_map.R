# table <- compound_mean_spagg |> dplyr::filter(class %in% c("Isoprene", "Monoterpenes", "Oxygenated-monoterpenes")) |> dplyr::select(-c(2,3,4))
# 
#table_pas_iso <- compound_mean_spagg |> dplyr::filter(class %in% c( "Monoterpenes", "Oxygenated-monoterpenes"))|> dplyr::select(-c(2,3,4))
# 
#  table_mono_emitters <- compound_mean_spagg |> dplyr::filter(class %in% c( "Monoterpenes", "Oxygenated-monoterpenes"))|> dplyr::select(-c(2,3,4)) |> dplyr::select(-c(
# "mean_fcar", "mean_acam" , "mean_msyl" ,  "mean_pang", "mean_pavi",
#  "mean_ppyr", "mean_ptre",  "mean_rala", "mean_rcat", "mean_spur", "mean_tbac", "mean_tgal",  "mean_scin",  "mean_ccog",
#  "mean_salb" , "mean_apse", "mean_pspi",  "mean_ainc", "mean_aopa" , "mean_cave" , "mean_faln" , "mean_iaqu",  "mean_pbri", "mean_pmah",  "mean_ppad" , "mean_sdom",
#  "mean_snig", "mean_stor", "mean_tafr", "mean_sari", "mean_srac",
# "mean_fexc", "mean_sauc", "mean_smou", "mean_sele"))
# 
# 
# 
# 
#  prefixes_a_supprimer <- c("fcar", "acam", "msyl", "pang", "pavi",
#                            "ppyr", "ptre", "rala", "rcat", "spur",
#                            "tbac", "tgal", "scin", "ccog", "salb",
#                            "apse", "pspi", "ainc", "aopa", "cave",
#                            "faln", "iaqu", "pbri", "pmah", "ppad",
#                            "sdom", "snig", "stor", "tafr", "sari",
#                            "srac", "fexc", "sauc", "smou", "sele")
# 
#  pattern <- paste0("^(", paste(prefixes_a_supprimer, collapse = "|"), ").*")
# 
#  table_mono_emitters <- compounds_samples_spagg_to_keep |>
#    dplyr::filter(class %in% c("Monoterpenes", "Oxygenated-monoterpenes")) |>
#    dplyr::select(-c(2, 3, 4)) |>
#    dplyr::select(-matches(pattern)) |>
#    dplyr::mutate(dplyr::across(-1, as.numeric))
# 
#  table_mono_emitters_chemodiv <- compounds_samples_spagg_to_keep |>
#    dplyr::filter(class %in% c("Isoprene", "Monoterpenes", "Oxygenated-monoterpenes")) |>
#    dplyr::select(-matches(pattern))
# 
# 
# 
# 
# 
# 
# 
# 
# 
# rename_mean_columns <- function(data, sp_screening) {
#   # Créer un vecteur de correspondance spcode -> nom scientifique avec espaces
#   name_mapping <- sp_screening |>
#     dplyr::select(spcode, full_scientific_name) |>
#     dplyr::distinct() |>
#     dplyr::mutate(
#       # Remplacer _ par espace
#       clean_name = gsub("_", " ", full_scientific_name),
#       # Nom de colonne original attendu (mean_ + spcode en minuscule)
#       old_col_name = paste0("mean_", tolower(spcode))
#     )
# 
#   # Créer un vecteur nommé pour le renommage
#   rename_vector <- setNames(name_mapping$clean_name, name_mapping$old_col_name)
# 
#   # Renommer les colonnes
#   for (i in seq_along(colnames(data))) {
#     col_name <- colnames(data)[i]
#     if (col_name %in% names(rename_vector)) {
#       colnames(data)[i] <- rename_vector[[col_name]]
#     }
#   }
# 
#   return(data)
# }
# 
# 
# 
# 
# library(ggplot2)
# library(dplyr)
# library(tidyr)
# library(scales)
# 
# # 1. Transformer les données
# table_long <- table |>
#   tidyr::pivot_longer(
#     cols = -compound,
#     names_to = "species",
#     values_to = "emission_rate"
#   ) |>
#   dplyr::mutate(
#     # Convertir en numérique
#     emission_rate_numeric = as.numeric(emission_rate),
#     # Remplacer NA par 0
#     emission_rate_numeric = ifelse(is.na(emission_rate_numeric), 0, emission_rate_numeric)
#   )
# 
# # 2. Calculer la proportion relative en pourcentage
# table_long <- table_long |>
#   dplyr::group_by(species) |>
#   dplyr::mutate(
#     total_species = sum(emission_rate_numeric, na.rm = TRUE),
#     relative_proportion = ifelse(total_species > 0,
#                                  (emission_rate_numeric / total_species) * 100,
#                                  0)
#   ) |>
#   ungroup()
# 
# # 3. Filtrer pour enlever les zéros (pas de bulle pour 0)
# table_filtered <- table_long |>
#   filter(emission_rate_numeric > 0)
# 
# # 4. Créer la bubble heatmap avec échelle log pour le taux d'émission
# ggplot(table_filtered,
#        aes(x = species,
#            y = compound,
#            size = log10(emission_rate_numeric + 1),  # +1 pour éviter -Inf quand emission=0
#            color = relative_proportion)) +
#   geom_point(alpha = 0.8) +
#   # Échelle de taille basée sur log10
#   scale_size_continuous(
#     name = "Taux d'émission (log10)",
#     range = c(1, 12),
#     breaks = c(0, 1, 2, 3),  # log10(1)=0, log10(10)=1, log10(100)=2, etc.
#     labels = c("1", "10", "100", "1000")  # Échelle inverse du log
#   ) +
#   # Échelle de couleurs pour le pourcentage
#   # scale_color_gradient(
#   #   name = "Proportion relative (%)",
#   #   colors = c(  "#ffb5b5", "#ff8181", "#db194a", "#c11239", "#8f061f"),
#   #   values = scales::rescale(c(0, 25, 50, 75, 100)),
#   #   limits = c(0, 100),
#   #   breaks = c(0, 25, 50, 75, 100)
#   # ) +
# 
#   scale_color_viridis_c(
#     name = "Proportion relative (%)",
#     option = "inferno",  # Palette inferno
#     direction = -1,       # -1 pour inverser les couleurs si besoin
#     limits = c(0, 100),
#     breaks = c(0, 25, 50, 75, 100)
#   ) +
# 
#   theme_minimal() +
#   theme(
#     axis.text.x = element_text(angle = 45, hjust = 1, size = 8),
#     axis.text.y = element_text(size = 7),
#     panel.grid.major = element_line(color = "gray90"),
#     legend.position = "right",
#     legend.box = "vertical"
#   ) +
#   labs(
#     x = "Espèces",
#     y = "Composés volatils",
#     title = "Carte de chaleur des émissions (échelle logarithmique)"
#   )
# 
