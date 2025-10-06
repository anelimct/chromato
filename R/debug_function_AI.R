trouver_fichiers_avec_cas <- function(files_path, numero_cas = "6753-98-6") {
  fichiers_trouves <- c()
  
  for (i in seq_along(files_path)) {
    data <- openxlsx::read.xlsx(files_path[i])
    
    # Nettoyage des données
    data <- data[-1, ]  
    colnames(data) <- data[1, ]  
    data <- data[-1, ]
    
    # Vérifier si le numéro CAS est présent
    if (numero_cas %in% data$`CAS Registry Number`) {
      fichiers_trouves <- c(fichiers_trouves, basename(files_path[i]))
    }
  }
  
  return(fichiers_trouves)
}

