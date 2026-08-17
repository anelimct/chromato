

compute_mean_EFtaxon_across_pop <- function(data, woodiv_species) {
  
  
  gragg_to_name <- woodiv_species |>
    dplyr::distinct(gragg, full_scientific_name) |>
    dplyr::group_by(gragg) |>
    # Prendre la première occurrence pour chaque gragg
    dplyr::slice(1) |>
    dplyr::ungroup()
  
  
  
  #Faire la moyenne de chaque population 
  
  pop_means <- data |>
    dplyr::group_by(gragg, Origin_pop, Compound) |>
    dplyr::summarise(
      EF_pop_mean = mean(EF, na.rm = TRUE),
      n_by_pop = dplyr::n(),  # Nombre d'observations par pop
      .groups = "drop"
    )
  ## Enregistrer l'objet R en enregistrer les deux sorties sous forme de liste puis selctionné la position del'objet souhaité dans les targets
  
  
  
  # Table avec isoprene et monoterpenes moyenné à l'échelle de la pop et pour les pop ou il y à la fois isoprene et monoterpenes compute aussi la sum à l'échelle de la pop  
  pop_wide <- pop_means |>
    dplyr::select(gragg, Origin_pop, Compound, EF_pop_mean) |> # enlever la colonne de n_by_pop
    tidyr::pivot_wider(names_from = Compound, values_from = EF_pop_mean) |>
    dplyr::mutate(
      sum_isoprenoids_pop = dplyr::if_else(
        !is.na(isoprene) & !is.na(monoterpenes),
        isoprene + monoterpenes,
        NA_real_
      )
    )
  
  
  
  # --- 2. Moyenne/min/max of sum à l'échelle espèce ---
  # (uniquement sur les pops où la somme a pu être calculée)
  species_sum <- pop_wide |>
    dplyr::filter(!is.na(sum_isoprenoids_pop)) |>
    dplyr::group_by(gragg) |>
    dplyr::summarise(
      sum_isoprenoids_mean = mean(sum_isoprenoids_pop, na.rm = TRUE),
      sum_isoprenoids_min  = min(sum_isoprenoids_pop, na.rm = TRUE),
      sum_isoprenoids_max  = max(sum_isoprenoids_pop, na.rm = TRUE),
      n_pop_sum = dplyr::n(),
      .groups = "drop"
    )
  
  
  # --- 3. species_means EF, moyenne des pop dispo iso et mono séparés = EF iso n'est pas forcément compute sur les mêmes pop que mono 
  # Moyenne par taxon across population = moyenne des moyennes de pop pour iso et pour mono
  # avec compte du nombre de pop utilisées
  species_means <- pop_means |>
    dplyr::group_by(gragg, Compound) |>
    dplyr::summarise(
      EF_species_mean = mean(EF_pop_mean, na.rm = TRUE),
      EF_species_min  = min(EF_pop_mean, na.rm = TRUE),
      EF_species_max  = max(EF_pop_mean, na.rm = TRUE),
      n_populations = dplyr::n(),
      .groups = "drop"
    )
  
  # --- 4. Mise en forme finale ---
  final_table <- species_means |>
    tidyr::pivot_wider(
      names_from = Compound,
      values_from = c(EF_species_mean, EF_species_min, EF_species_max, n_populations),
      names_glue = "{Compound}_{.value}"
    ) |>
    dplyr::rename(
      isoprene = isoprene_EF_species_mean,
      isoprene_min = isoprene_EF_species_min,
      isoprene_max = isoprene_EF_species_max,
      monoterpenes = monoterpenes_EF_species_mean,
      monoterpenes_min = monoterpenes_EF_species_min,
      monoterpenes_max = monoterpenes_EF_species_max,
      n_pop_isoprene = isoprene_n_populations,
      n_pop_monoterpenes = monoterpenes_n_populations
    ) |>
    dplyr::left_join(species_sum, by = "gragg") |>          # <- ajout de la somme
    dplyr::left_join(gragg_to_name, by = "gragg") |>
    dplyr::select(
      "full_scientific_name", "gragg",
      "isoprene", "isoprene_min", "isoprene_max",
      "monoterpenes", "monoterpenes_min", "monoterpenes_max",
      "sum_isoprenoids_mean", "sum_isoprenoids_min", "sum_isoprenoids_max",
      "n_pop_isoprene", "n_pop_monoterpenes", "n_pop_sum"
    ) |>
    dplyr::rename("name_complete" = "full_scientific_name") |>
    dplyr::mutate(name_complete = dplyr::case_when(
      name_complete == "Juniperus_deltoides" ~ "Juniperus_oxycedrus",
      TRUE ~ name_complete
    ))
  
  return(final_table)
}


# compute_DB_bvocs_iso_mono_EF <- function(all_data_mean_EF_taxon, type){
#   
#   if(type == "figure") {
#     DB_bvocs_iso_mono_EF <- all_data_mean_EF_taxon |> 
#       tibble::column_to_rownames(var = "name_complete") |>
#       dplyr::mutate(Sum = isoprene + monoterpenes) |>
#       dplyr::select("isoprene", "monoterpenes", "Sum")
#   }
#   else {
#     DB_bvocs_iso_mono_EF <- all_data_mean_EF_taxon |> 
#       tibble::column_to_rownames(var = "gragg") |>
#       dplyr::mutate(Sum = isoprene + monoterpenes) |>
#       dplyr::select("isoprene", "monoterpenes", "Sum")
#   }
#   
#   return(DB_bvocs_iso_mono_EF)
# }

normaliser_dataframe <- function(dataframe) {
  # Calculer les moyennes et les écarts-types de chaque variable
  means <- apply(dataframe, MARGIN = 2, FUN = function(x) mean(x, na.rm = TRUE))
  sds <- apply(dataframe, MARGIN = 2, FUN = function(x) sd(x, na.rm = TRUE))
  
  # Normaliser les données en utilisant les moyennes et les écarts-types calculés
  dataframe_normalise <- scale(dataframe, center = means, scale = sds) |>  
    as.data.frame()
  
  return(dataframe_normalise)
}


merge_trait_EF <- function(imputed.traits_3T, DB_bvocs_iso_mono_EF){
  
  merged_data_3T <- merge(imputed.traits_3T, DB_bvocs_iso_mono_EF, by = "row.names", all.x = TRUE) |>
    tibble::column_to_rownames(var = "Row.names") |>
    dplyr::mutate(
      total = sum_isoprenoids_mean,          # somme calculée au niveau pop, plus fiable que isoprene+monoterpenes recalculé ici
      p_isoprene = isoprene / total,
      p_monoterpenes = monoterpenes / total,
      prct_isoprene = p_isoprene * 100,
      prct_monoterpenes = p_monoterpenes * 100
    ) |>
    dplyr::mutate(isoprene_mod = ifelse(isoprene == 0, 0.00001, isoprene)) |>
    dplyr::mutate(BVOCsData = dplyr::case_when(
      is.na(isoprene) & is.na(monoterpenes) ~ "0",   # aucune mesure BVOC du tout
      .default = "1"
    )) |>
    dplyr::mutate(type = dplyr::case_when(
      is.na(isoprene) & is.na(monoterpenes) ~ NA_character_,
      isoprene >= 1 & monoterpenes > 0.2 ~ "both",
      monoterpenes > 0.2 ~ "mono",
      isoprene > 1 ~ "iso",
      .default = "NE"
    )) |>
    dplyr::mutate(binaire = dplyr::case_when(
      type == "iso" ~ 1,
      type == "mono" ~ 0,
      .default = NA_real_
    ))
  
  return(merged_data_3T)
}
  
create_residual_correlogram <- function(tree, residuals, col_name = "Residuals") {
  # Vérifier le type de données des résidus
  if (!is.vector(residuals)) {
    stop("Les résidus doivent être un vecteur.")
  }
  
  # Créer un objet phylo4 avec l'arbre élagué
  pruned_tree <- ape::drop.tip(tree, setdiff(tree$tip.label, names(residuals)))
  phylo4_object <- phylobase::phylo4(pruned_tree)
  
  # Créer un data.frame avec les résidus et nommer la colonne
  df_residuals <- data.frame(residuals)
  colnames(df_residuals) <- col_name
  
  # Créer un objet phylo4d avec les résidus
  phylo4d_object <- phylobase::phylo4d(phylo4_object, tip.data = df_residuals)
  
  # Créer le correlogramme
  correlogram <- phylosignal::phyloCorrelogram(
    phylo4d_object,
    trait = col_name,
    dist.phylo = "patristic"
  )
  
  mesures <- phylosignal::phyloSignal(phylo4d_object, methods = "all")
  print(mesures)
  
  # Afficher le plot
  phylosignal:: plot.phylocorrelogram(correlogram, main = paste("Phylogenetic correlogram of", col_name))
  
  return(correlogram)
}



compute_total_ITV <- function(data, woodiv_species) {
  
  gragg_to_name <- woodiv_species |>
    dplyr::distinct(gragg, full_scientific_name) |>
    dplyr::group_by(gragg) |>
    dplyr::slice(1) |>
    dplyr::ungroup()
  
  # Total variance across all individuals (all populations pooled)
  species_stats <- data |>
    dplyr::group_by(gragg, Compound) |>
    dplyr::summarise(
      EF_mean = mean(EF, na.rm = TRUE),
      EF_var_total = var(EF, na.rm = TRUE),   # <-- total variance
      EF_sd_total = sd(EF, na.rm = TRUE),     # standard deviation for interpretability
      n_total = dplyr::n(),                    # total number of measurements
      .groups = "drop"
    )
  
  # Pivot to wide format (isoprene / monoterpenes)
  final_table <- species_stats |>
    tidyr::pivot_wider(
      names_from = Compound,
      values_from = c(EF_mean, EF_var_total, EF_sd_total, n_total),
      names_glue = "{Compound}_{.value}"
    ) |>
    dplyr::rename(
      isoprene_mean = `isoprene_EF_mean`,
      monoterpenes_mean = `monoterpenes_EF_mean`,
      isoprene_var = `isoprene_EF_var_total`,
      monoterpenes_var = `monoterpenes_EF_var_total`,
      isoprene_sd = `isoprene_EF_sd_total`,
      monoterpenes_sd = `monoterpenes_EF_sd_total`,
      n_isoprene = `isoprene_n_total`,
      n_monoterpenes = `monoterpenes_n_total`
    ) |>
    dplyr::left_join(gragg_to_name, by = "gragg") |>
    dplyr::select(
      name_complete = full_scientific_name, 
      gragg, 
      isoprene_mean, monoterpenes_mean,
      isoprene_var, monoterpenes_var,
      isoprene_sd, monoterpenes_sd,
      n_isoprene, n_monoterpenes
    ) |>
    dplyr::mutate(name_complete = dplyr::case_when(
      name_complete == "Juniperus_deltoides" ~ "Juniperus_oxycedrus",
      TRUE ~ name_complete
    ))
  
  return(final_table)
}