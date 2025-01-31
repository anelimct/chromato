#' See if files present in the chromato file for the year match files said to be present
#'
#' @param data A dataframe with bvocs samples 
#' @param year 
#'
#' @return A dataframe with : how many chromato should by present, how many are available, the names the missing files if there is any within pairs of samples ( for the year 2023 one tube C5 one tube C15)
#' @export
#'
#' @examples
chromato_file_check <- function(data, year) {
  # # Identifier les paires de tubes
  data <- data |>
    dplyr::group_by(Date, Taxon) |>
    dplyr::mutate(pair_id = dplyr::cur_group_id()) |>
    dplyr::ungroup()
  # 
  # # Liste des fichiers dans le dossier
  files_in_folder <- list.files(here::here("data", "GC-MS_files", paste0(year)))
  # 
  # # Ajouter une colonne available
  data <- data |>
    dplyr::mutate(available = ifelse(paste0(ID, ".D") %in% files_in_folder, TRUE, FALSE))|>
    dplyr::mutate(available = as.logical(available))
  
  data <- data |>
    dplyr::group_by(pair_id, Taxon, Date)|>
    dplyr::summarize( theorique_count = dplyr::n(), available_count = sum(available == "TRUE"), compo_list = paste0(Composés[available == "TRUE"], collapse = ", "), present_list = paste0(ID[available == "TRUE"], collapse = ", "), missing_list = paste0(ID[available == "FALSE"], collapse = ", "))|>
    dplyr::ungroup()
  
  return(data)
}



#' Creates folders with chromato organised by batch (same week of desorption)
#'
#' @param data 
#' @param year 
#'
#' @return
#' @export
#'
#' @examples
organize_gc_files_by_batch <- function(data, year) {
  
  
blanks <- data |> 
      dplyr::select(batch,`list(ID_blank)` ) |> 
       dplyr::filter(!is.na(batch), !purrr::map_lgl(`list(ID_blank)`, is.null)) |>  # Filtrer les lignes pertinentes
       tidyr::unnest_longer(`list(ID_blank)`) |>                            # Éclater les éléments des listes
    dplyr::rename(ID = `list(ID_blank)`) |> dplyr::distinct(ID, batch)

data <- data |>  dplyr::filter(!startsWith(ID, "B_")) |>  dplyr::select(ID, batch) |>  rbind(blanks)

  # Définir les dossiers source et cible
  source_folder <- here::here("data", "GC-MS_files", paste0(year))
  target_base_folder <- here::here("data", "GC-MS_batches", paste0(year))
  
  # Liste des fichiers/dossiers dans le dossier source
  files_in_folder <- list.files(source_folder, full.names = TRUE)
  
  # Ajouter une colonne 'available' pour vérifier l'existence des fichiers/dossiers
  data <- data |>
    dplyr::mutate(
      available = paste0(ID, ".D") %in% basename(files_in_folder)
    )
  
  # Filtrer les lignes valides
  data <- data |>
    dplyr::filter(!is.na(batch) & available)
  
  # Ajouter les chemins source et cible
  data <- data |>
    dplyr::mutate(
      source_file = file.path(source_folder, paste0(ID, ".D")),
      target_folder = file.path(target_base_folder, as.character(batch))
    )
  
  # Vérifier ou créer les dossiers de batch
  data <- data |> 
    dplyr::rowwise() |>
    dplyr::mutate(
      folder_created = {
        if (!dir.exists(target_folder)) {
          dir.create(target_folder, recursive = TRUE)
        }
        target_folder
      }
    ) |>
    dplyr::ungroup()
  
  # Copier les fichiers/répertoires
  data <- data |>
    dplyr::mutate(
      copy_success = purrr::map2_lgl(source_file, target_folder, ~ {
        if (file.exists(.x)) {
          if (dir.exists(.x)) { # Cas où la source est un répertoire
            file.copy(.x, .y, overwrite = TRUE, recursive = TRUE)
          } else { # Cas où la source est un fichier
            file.copy(.x, file.path(.y, basename(.x)), overwrite = TRUE)
          }
        } else {
          FALSE
        }
      })
    )
  
  # Afficher les fichiers ou répertoires non copiés
  failed_copies <- data |> dplyr::filter(!copy_success)
  if (nrow(failed_copies) > 0) {
    warning("Certains fichiers ou répertoires n'ont pas été copiés :\n", 
            paste(failed_copies$source_file, collapse = "\n"))
  }
  
  # Retourner les données mises à jour
  return(data)
}




#' Rename alcanes files in the data folder
#'
#' @param alcanes 
#'
#' @return
#' @export
#'
#' @examples
rename_GC_files <- function(dataframe, name){
  source_folder <- here::here("data", "GC-MS_files", paste0(name))
  files_in_folder <- list.files(source_folder, full.names = FALSE)
  
  # Parcourir les correspondances et renommer les fichiers/dossiers
  for (i in 1:nrow(dataframe)) {
    initial <-dataframe$initial_name[i]
    new <- dataframe$new_name[i]
    # Chemins source et cible
    initial_path <- file.path(source_folder, initial)
    new_path <- file.path(source_folder, new)
    
    # Renommer le dossier
    if (file.exists(initial_path)) {
      success <- file.rename(initial_path, new_path)
      if (!success) {
        warning(paste("Impossible de renommer :", initial, "->", new))
      }
    } else {
      warning(paste("Le fichier/dossier", initial, "n'existe pas dans le dossier source."))
    }
  }
  return(dataframe)
}


