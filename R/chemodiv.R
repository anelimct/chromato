format_sample_data_chemodiv <- function(compounds_table, times_compound_sp, valid_samples_mono){
  
  valid_samples <- stringr::str_to_lower(valid_samples_mono$ID)
  valid_samples <- stringr::str_c(valid_samples, ".cdf") 
  
  
  first<- as.data.frame(compounds_table)
  
  first_transposed <- as.data.frame(t(first))
  
  # Renommer les colonnes avec les compounds (première ligne du tableau original)
  colnames(first_transposed) <- first[, 1]
  
  # Supprimer la première ligne qui contient maintenant les noms des compounds
  first_transposed <- first_transposed[-1, ]
  
  # Convertir toutes les colonnes en numeric
  first_transposed[] <- lapply(first_transposed, function(x) as.numeric(as.character(x)))
  
  # Supprimer les colonnes qui n'ont que des NA
  first_transposed <- first_transposed[, colSums(is.na(first_transposed)) < nrow(first_transposed)]
  
  # Ajouter les noms des samples comme colonne
  first_transposed$sample <- rownames(first_transposed)
  
  # Réorganiser pour avoir sample en première colonne
  first_transposed <- first_transposed[, c("sample", setdiff(colnames(first_transposed), "sample"))]
  
  
  
  
  sampledata <- first_transposed[, -1] 
  sampledata <- sampledata[rowSums(!is.na(sampledata)) > 0, ]
  sampledata[is.na(sampledata)] <- 0 
  
  
  
  
  ## 
  
  compound_data <- compounds_table[, c(1:4)]
  

  
  # Get the current rownames of sampledata
  current_rownames <- rownames(sampledata)
  
  # Find which rownames are in valid_samples
  matching_rows <- current_rownames %in% valid_samples
  
  # Keep only the rows where rownames match valid_samples
  sampledata <- sampledata[matching_rows, ]
  
  ## filtre pour ne garder des selement les samples des espèces qui ont étét échantilonné au moind  trois fois 

  
  
  spagg <- unique(times_compound_sp$spagg)
  sampledata <- sampledata[substr(rownames(sampledata), 1, 4) %in% spagg, ]
  ##Il faudrait aussi enlever tous les componds qui sont des sesquiterpenes comme il n'ont pa spu être quantifié pour tous 
  comp_to_keep <- compound_data |>  dplyr::filter(class %in% c("Isoprene", "Monoterpenes", "Oxygenated-monoterpenes" ) )
  
  comp_to_keep <- unique(comp_to_keep$compound)
  
  compound_data <- compound_data |>  dplyr::filter(compound %in% comp_to_keep)
  sampledata <- sampledata[, colnames(sampledata) %in% comp_to_keep]
  
  sampledata_relative <- as.data.frame(t(apply(sampledata, 1, function(x) x/sum(x))))
  
  ##
  
  
  pheno <- WOODIV_v2_Trait_data |> 
    dplyr::filter(trait == "LeafPheno") |> 
    dplyr::rename("spagg" = "spcode")
  
  group_df <- data.frame(
    sample = rownames(sampledata_relative),
    species = substr(rownames(sampledata_relative), 1, 4)
  ) |>  dplyr::mutate(spagg = stringr::str_to_upper(species)) |> dplyr::left_join( pheno)
  
  
  group <- group_df[,3]
  
  
  
}
