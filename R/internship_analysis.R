

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
      n_pop = dplyr::n(),  # Nombre d'observations par pop
      .groups = "drop"
    )
  
  # Moyenne par taxon across population = moyenne des moyennes de pop
  # avec compte du nombre de pop utilisées
  species_means <- pop_means |>
    dplyr::group_by(gragg, Compound) |>
    dplyr::summarise(
      EF_species_mean = mean(EF_pop_mean, na.rm = TRUE),
      n_populations = dplyr::n(),  # Nombre de populations pour ce Taxon-Compound
      .groups = "drop"
    )
  
  #Mise en forme du tableau
  final_table <- species_means |>
    tidyr::pivot_wider(
      names_from = Compound,
      values_from = c(EF_species_mean, n_populations),
      names_glue = "{Compound}_{.value}"
    ) |>
    # Renommer les colonnes pour plus de clarté
    dplyr::rename(
      isoprene = `isoprene_EF_species_mean`,
      monoterpenes = `monoterpenes_EF_species_mean`,
      n_pop_isoprene = `isoprene_n_populations`,
      n_pop_monoterpenes = `monoterpenes_n_populations`
    ) |> dplyr::left_join(  gragg_to_name, by = "gragg") |> dplyr::select("full_scientific_name", "gragg", "isoprene", "monoterpenes","n_pop_isoprene", "n_pop_monoterpenes") |> dplyr::rename("name_complete"= "full_scientific_name" )
  
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
  dataframe_normalise <- scale(dataframe, center = means, scale = sds) %>% 
    as.data.frame()
  
  return(dataframe_normalise)
}


merge_trait_EF<- function(imputed.traits_3T, DB_bvocs_iso_mono_EF){
  
  merged_data_3T <- merge(imputed.traits_3T, DB_bvocs_iso_mono_EF, by = "row.names", all.x = TRUE) |>  
    tibble::column_to_rownames(var = "Row.names") |> 
   # merge(splist_storing_2_, by ="row.names", all.x = TRUE) |> 
    #tibble::column_to_rownames(var = "Row.names") |> 
    dplyr::mutate(total = isoprene + monoterpenes,
           p_isoprene = (isoprene / total),
           p_monoterpenes = (monoterpenes / total), 
           prct_isoprene = (isoprene / total) * 100, 
           prct_monoterpenes = (monoterpenes / total) * 100) |> 
    dplyr::mutate(isoprene_mod = ifelse(isoprene == 0, 0.00001, isoprene)) |> 
    dplyr::mutate(BVOCsData = dplyr::case_when(is.na(Sum) == TRUE ~ "0",
                                 .default = "1")) |> 
    dplyr::mutate(type = dplyr::if_else(!is.na(Sum), 
                          ifelse(isoprene >= 1 & monoterpenes > 0.2, "both", 
                                 ifelse(monoterpenes > 0.2, "mono", 
                                        ifelse(isoprene > 1, "iso", "NE"))), 
                          NA_character_)) |> 
    dplyr::mutate(binaire = ifelse(type == "iso", 1, ifelse(type == "mono", 0, NA))) 
    # dplyr::mutate(nouveau_type = dplyr::if_else(type == "iso", type,
    #                               dplyr::if_else(type == "NE", type,
    #                                       dplyr::if_else(Stockage == "non", type,
    #                                               dplyr::if_else(Stockage == "oui", paste0(type, "_s"), NA_character_)))))
}

  
create_residual_correlogram <- function(tree, residuals, col_name = "Residuals") {
  # Vérifier le type de données des résidus
  if (!is.vector(residuals)) {
    stop("Les résidus doivent être un vecteur.")
  }
  
  # Créer un objet phylo4 avec l'arbre élagué
  pruned_tree <- ape::drop.tip(tree, setdiff(tree$tip.label, names(residuals)))
  phylo4_object <- as(pruned_tree, "phylo4")
  
  # Créer un data.frame avec les résidus et nommer la colonne
  df_residuals <- data.frame(residuals)
  colnames(df_residuals) <- col_name
  
  # Créer un objet phylo4d avec les résidus
  phylo4d_object <- phylo4d(phylo4_object, tip.data = df_residuals)
  
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
  
}
