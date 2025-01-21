create_library <- function(files_path){
  CAS_lib <- c()
  
  for (i in seq_along(files_path)) {
    data <- openxlsx::read.xlsx(files_path[i])
    
    
    data <- data[-1, ]  
    colnames(data) <- data[1, ]  #
    data <- data[-1, ]
    
    CAS_lib <- rbind(CAS_lib, data)
    CAS_lib <- unique(CAS_lib)
    
  }
  return(CAS_lib)
}
