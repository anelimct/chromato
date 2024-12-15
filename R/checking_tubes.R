#' See if files present in the chromato file for the year match files said to be present
#'
#' @param data 
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


organize_gc_files_by_batch <- function(data, year) {
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