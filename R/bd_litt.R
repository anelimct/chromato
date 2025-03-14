
read_excel_articles <- function(file) {
  readxl::read_xlsx(file, sheet = 2) 
}


#WOODIV_Species_code <- WOODIV_Species_code |> mutate(genus_species_subspecies = gsub('_', ' ', genus_species_subspecies))
#data <- WOODIV_Species_code %>%   left_join(WOODIV_Nomenclature %>% select(spcode, class, subclass, order, family), by = 'spcode')
#data <- data %>%    left_join(WOODIV_Aggregation , by = 'spcode')

species_aggregation <- function(data, species_woodiv){
  
  data <- data |> 
    dplyr::left_join(species_woodiv %>% dplyr::select(genus_species_subspecies, to_aggregate_with),
              by = c('Taxon_in_ref' = 'genus_species_subspecies')) |> 
    dplyr::rename(spcode.agg = to_aggregate_with)
}


select_iso_mono <- function(data){
  variantes_monoterpenes <- c("monoterpene", "monoterpens", "monterpene", "monterpenes", "monoterp", "monoterpeness")
  
  # Standardize Compound names and filter relevant compounds
  data <- data |> 
    dplyr::mutate(Compound = ifelse(Compound %in% variantes_monoterpenes, "monoterpenes", Compound)) |> 
    dplyr::filter(Compound == "isoprene" | Compound == "monoterpenes")
}

numeric_emissions_g_h <- function(data){
  
  data <- data |> 
  #traces and nd as 0
    dplyr::mutate(Emission = ifelse(Emission== "nd", "0", Emission)) |> 
    dplyr::mutate(Emission = ifelse(Emission== "tr", "0", Emission)) |>  #a voir avec Caroline
  #numeric emissions
    
    dplyr::mutate(Emission = as.numeric(Emission), 
           Emission_var_value = as.numeric(Emission_var_value)) |> 
    #emisions in hours
    dplyr::mutate(Emission = ifelse(Emission_unit_time == "s", Emission * 3600, Emission),
           Emission_var_value = ifelse(Emission_unit_time == "s", Emission_var_value * 3600, Emission_var_value),
           Emission_unit_time = ifelse(Emission_unit_time == "s", "h", Emission_unit_time)) |> 
    
    dplyr::mutate(Emission = ifelse(Emission_unit_comp == "nanog", Emission / 1000, Emission),
           Emission_var_value = ifelse(Emission_unit_comp == "nanog", Emission /1000, Emission_var_value),
           Emission_unit_comp = ifelse(Emission_unit_comp == "nanog", "microg", Emission_unit_comp)) |> 
    
    dplyr::mutate(Emission = ifelse(Emission_unit_comp == "g", Emission *1000, Emission),
           Emission_var_value = ifelse(Emission_unit_comp == "g", Emission *1000, Emission_var_value),
           Emission_unit_comp = ifelse(Emission_unit_comp == "g", "microg", Emission_unit_comp)) |> 
    
    dplyr::mutate(Emission = ifelse(Emission_unit_comp == "nanogC", Emission / 1000, Emission),
           Emission_var_value = ifelse(Emission_unit_comp == "nanogC", Emission /1000, Emission_var_value),
           Emission_unit_comp = ifelse(Emission_unit_comp == "nanogC", "microgC", Emission_unit_comp))
}


convert_temperature <- function(data) {
  # Si la température est en Celsius
  data <- data |> 
    dplyr::mutate(
           Temperature_C = dplyr::if_else(Temperature_unit == "K", Temperature - 273.15, Temperature),
           Temperature_min_C = dplyr::if_else(Temperature_unit == "K", Temperature_min - 273.15, Temperature_min),
           Temperature_max_C = dplyr::if_else(Temperature_unit == "K", Temperature_max - 273.15, Temperature_max))
  
  return(data)
}



select_std_or_standardisable <- function(data, sp_storing) {
  
  data <- data |> 
    dplyr::left_join(sp_storing, by = c('spcode.agg' = 'spcode' ))
  
  data_std <- data |> dplyr::filter(Standardized == "true")
  
  data<- data |> dplyr::filter(Standardized == "false")
  
  data <- data |> dplyr::filter(
    Stockage == "oui" & Compound == "monoterpenes" & 
      (Temperature != "NA" | (Temperature_min != "NA" & Temperature_max != "NA")) |
      
    Stockage == "non"& Compound == "monoterpenes" &  ((Temperature != "NA" & PAR != "NA") | 
      
      (Temperature_min != "NA" & Temperature_max != "NA" & PAR_min != "NA" & PAR_max != "NA")  | 
      (Temperature != "NA"& PAR_min != "NA" & PAR_max != "NA") |
      (PAR != "NA" & Temperature_min != "NA" & Temperature_max != "NA")) |
      
      Compound == "isoprene" & 
      ((Temperature != "NA" & PAR != "NA") | 
      (Temperature_min != "NA" & Temperature_max != "NA" & PAR_min != "NA" & PAR_max != "NA")  | 
      (Temperature != "NA"& PAR_min != "NA" & PAR_max != "NA") |
      (PAR != "NA" & Temperature_min != "NA" & Temperature_max != "NA"))) 
  
  data_combined <- rbind(data_std, data)
    
} 


select_months <- function(data, start_month, end_month) {
  # Convert arguments to numeric
  start_month_num <- as.numeric(start_month)
  end_month_num <- as.numeric(end_month)
  
  # Convert mm_1 and mm_2 to numeric, handling NA values
  data <- data |> 
    dplyr::mutate(
      mm_1 = as.numeric(mm_1),
      mm_2 = as.numeric(mm_2)
    )
  
  # Filter the data
  data <- data |> 
    dplyr::filter(
      (mm_1 >= start_month_num & mm_1 <= end_month_num) |
        (mm_2 >= start_month_num & mm_2 <= end_month_num) |
        (is.na(mm_1) & is.na(mm_2))
    )
  
  return(data)
}
 




select_temp_and_par <- function(data, Tmax, Tmin, PARmax, PARmin){
  
  
  
  
  data <-   data |>  dplyr::mutate(dplyr::across(c(PAR, PAR_min, PAR_max, Temperature, Temperature_min, Temperature_max), as.numeric))  |>
    convert_temperature() |> 
    dplyr::filter(
      (is.na(Temperature_C) | (Temperature_C >= Tmin & Temperature_C <= Tmax)) &
        (is.na(Temperature_min_C) | Temperature_min_C >= Tmin) &
        (is.na(Temperature_max_C) | Temperature_max_C <= Tmax)
    ) |>  
    dplyr::filter(
      (is.na(PAR) | (PAR >= PARmin & PAR <= PARmax)) &
        (is.na(PAR_min) | PAR_min >= PARmin) &
        (is.na(PAR_max) | PAR_max <= PARmax)
    )
    #mutate(Temperature_range = Temperature_max - Temperature_min,PAR_range = PAR_max - PAR_min)
}




create_venn_diagram <- function(data) {
  # Subset the data for monoterpenes and isoprene
  monoterpenes <- dplyr::filter(data, Compound %in% c("monoterpenes"))
  isoprene <- dplyr::filter(data, Compound %in% c("isoprene"))
  
  # Remove NA values from spcode
  monoterpenes <- na.omit(monoterpenes$spcode.agg)
  isoprene <- na.omit(isoprene$spcode.agg)
  
  # Get unique species for each category
  species_monoterpenes <- unique(monoterpenes)
  species_isoprene <- unique(isoprene)
  
  # Create the Venn diagram
  venn.plot <- VennDiagram::venn.diagram(
    x = list(Monoterpenes = species_monoterpenes, Isoprene = species_isoprene),
    category.names = c("Monoterpenes", "Isoprene"),
    filename = NULL,
    output = TRUE,
    height = 480,
    width = 480,
    resolution = 300,
    col = c("skyblue", "pink"),
    fill = c(scales::alpha("skyblue", 0.5), scales::alpha("pink", 0.5)),
    fontfamily = "sans",
    cat.cex = 0.8,
    cat.fontfamily = "sans",
    cat.col = c("black", "black"),
    cat.dist = c(0.07, 0.07),
    margin = 0.05
  )
  
  # Draw the Venn diagram
  grid::grid.draw(venn.plot)
}


#### Population Variable---



clean_country <- function(country) {
  country <- ifelse(country != "NA", stringr::str_to_title(trimws(country)), country)
  country <- ifelse(country == "Iltaly", "Italy", country)
  return(country)
}

# Function to remove accents
remove_accents <- function(text) {
  text <- chartr("àâäéèêëîïôöùûüç", "aaaeeeeiioouuuc", text)
  return(text)
}

# Function to clean and format city and locality names
clean_city_locality <- function(name) {
  name <- remove_accents(name)
  name <- ifelse(name != "NA", stringr::str_to_title(trimws(name)), name)
  name <- gsub("\\b(Castelporzi(a|ano)\\b|Castelporziano Nature Reserve)", "Castelporziano", name)
  return(name)
}



create_population_variable <- function(df) {
  # Liste des groupes similaires
  groupes_similaires <- list(
    "Burrina" = c("Burrina Cirat", "Burrina NA", "Cirat NA"),
    "Tartu" = c("Tartu NA", "NA Tartu"),
    "Castelporziano" = c("Castelporziano NA", "Castelporziano Castelporziano", "Castelporziano Santo Querico"),
    "Freiburg" = c("Freiburg Hartheimer Wald", "NA Hartheimer, Near Freiburg"),
    "Madrid" = c("Madrid Tres Cantos", "Madrid NA"),
    "Julich" = c("Julich Julich Research Center", "Julich NA"),
    "Hyytiala" = c("NA Hyytiala", "Hyytiala Station For Measuring Forest Ecosystem–Atmosphere Relations", "Hyytiala NA"),
    "Arad" = c("Arad Lipova Forest", "Arad NA")
  )
  
  
  # Create the group dataframe
  df_groupes <- data.frame(
    Origin_pop = unlist(groupes_similaires),
    nom_groupe = rep(names(groupes_similaires), sapply(groupes_similaires, length))
  )
  
  # Merge the dataframes
  df <- dplyr::left_join(df, df_groupes, by = "Origin_pop")
  
  # Update the Origin_pop column
  df$Origin_pop <- ifelse(!is.na(df$nom_groupe), df$nom_groupe, df$Origin_pop)
  
  # Remove the nom_groupe column
  df <- df |>   dplyr::select(-nom_groupe)
  
  return(df)
}



