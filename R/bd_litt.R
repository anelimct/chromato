
read_excel_articles <- function(file) {
  readxl::read_xlsx(file, sheet = 1) 
}


#WOODIV_Species_code <- WOODIV_Species_code |> mutate(genus_species_subspecies = gsub('_', ' ', genus_species_subspecies))
#data <- WOODIV_Species_code %>%   left_join(WOODIV_Nomenclature %>% select(spcode, class, subclass, order, family), by = 'spcode')
#data <- data %>%    left_join(WOODIV_Aggregation , by = 'spcode')

species_aggregation <- function(data, species_woodiv, dataset){
  if (dataset == "litt"){
    
    species_woodiv <- species_woodiv |> 
      dplyr::group_by(gragg) |> 
      dplyr::mutate(
        taxon = dplyr::first(full_scientific_name),
        Taxon = stringr::str_replace_all(taxon, "_", " ")
      ) |> 
      dplyr::ungroup() 
    
    #problematic subspecies that fuck up the joint
    name_mapping <- tibble::tribble(
      ~Taxon_in_ref, ~corrected_name,
      "Quercus macrolepis", "Quercus ithaburensis macrolepis",
      "Quercus ilex subsp. ilex",  "Quercus ilex ilex",
      "Quercus ilex subsp. rotundifolia", "Quercus ilex rotundifolia", 
      "Quercus rotundifolia", "Quercus ilex rotundifolia"
      # Add other problematic mappings here
    )
    
    
    data <- data |> dplyr::filter(!is.na(Taxon_in_ref )) |> 
      dplyr::left_join(name_mapping, by = "Taxon_in_ref") |>
      dplyr::mutate(
        match_name = dplyr::coalesce(corrected_name, Taxon_in_ref)
      ) |>
      dplyr::left_join(
        species_woodiv |>
          dplyr::select(full_scientific_name, spagg, gragg, Taxon) |>
          dplyr::mutate(full_scientific_name = stringr::str_replace_all(full_scientific_name, "_", " ")),
        by = c('match_name' = 'full_scientific_name')
      ) |>
      dplyr::select(-corrected_name, -match_name) |> dplyr::mutate(
        Taxon = stringr::str_replace_all(Taxon, " ", "_"))
    
  } 
  else if (dataset == "field_"){
    data <- data |> 
      dplyr::left_join(species_woodiv |>  dplyr::select(full_scientific_name, spagg, gragg),
                       by = c('Taxon' = 'full_scientific_name'))
    
  }
  
  else {
    data <- data |> 
      dplyr::left_join(species_woodiv |>  dplyr::select(full_scientific_name, spagg, gragg) |> dplyr::mutate( full_scientific_name = stringr::str_replace_all(full_scientific_name, "_", " ") ),
                       by = c('Taxon' = 'full_scientific_name'))
  }
  
}


select_iso_mono <- function(data){
  variantes_monoterpenes <- c("monoterpene", "monoterpens", "monterpene", "monterpenes", "monoterp", "monoterpeness")
  
  # Standardize Compound names and filter relevant compounds
  data <- data |> 
    dplyr::mutate(Compound = ifelse(Compound %in% variantes_monoterpenes, "monoterpenes", Compound)) |> 
    dplyr::filter(Compound == "isoprene" | Compound == "monoterpenes"  | Compound == "sesquiterpenes" | Compound == "ox-sesquiterpenes"| Compound == "ox-monoterpenes") |> 
    dplyr::filter(Compound_TvsP != "partial")
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
  
  
  data <- data |>  dplyr:: filter(Emission_unit_leaf == "g")
}

compound_unit_µg <- function(data) {
  ## passage de microgC à microC, pour cela il faut la fraction moléaire du carbone dans la molécule
  # pour monoterpenes C10H16 = 12.01*10 + 16 = 136.1 g.mol, fraction molaire carbone = 120.1 , fraction massique  = 0.8824394, Emission µgC/0.8824394 = Emission µg
  # pour isoprene C5H8 = 12.01*5 + 8 = 68.05 g.mol, fraction molaire carbone = 60.05 , fraction massique  = 0.8824394, Emission µgC/0.8824394 = Emission µg
  data <- data |>   dplyr::mutate(Emission = ifelse(Emission_unit_comp == "microgC" , Emission/0.8824394 , Emission)) 
  data <- data |>   dplyr::mutate(Emission = ifelse(Emission_unit_comp == "nanomol" & Compound == "monoterpenes", Emission*136.1*1e+06, Emission ))
  data <- data |>   dplyr::mutate(Emission = ifelse(Emission_unit_comp == "nanomol" & Compound == "isoprene", Emission*68.05*1e+06, Emission ))
  
  
  data <- data |>   dplyr::mutate(Emission_var_value = ifelse(Emission_unit_comp == "microgC" , Emission_var_value/0.8824394 , Emission_var_value)) 
  data <- data |>   dplyr::mutate(Emission_var_value = ifelse(Emission_unit_comp == "nanomol" & Compound == "monoterpenes", Emission_var_value*136.1*1e+06, Emission_var_value ))
  data <- data |>   dplyr::mutate(Emission_var_value = ifelse(Emission_unit_comp == "nanomol" & Compound == "isoprene", Emission_var_value*68.05*1e+06, Emission_var_value ))
  
  data$Emission_unit_comp <- "microg"
  return(data)
}



convert_temperature <- function(data) {
  # Si la température est en Celsius
  data <- data |> 
    dplyr::mutate(Temperature_C = dplyr::if_else(Temperature_unit == "C", Temperature, Temperature - 273.15),
           Temperature_min_C = dplyr::if_else(Temperature_unit == "C", Temperature_min, Temperature_min - 273.15),
           Temperature_max_C = dplyr::if_else(Temperature_unit == "C", Temperature_max, Temperature_max - 273.15),
           Temperature_K = dplyr::if_else(Temperature_unit == "C", Temperature + 273.15, Temperature),
           Temperature_min_K = dplyr::if_else(Temperature_unit == "C", Temperature_min + 273.15, Temperature_min),
           Temperature_max_K = dplyr::if_else(Temperature_unit == "C", Temperature_max + 273.15, Temperature_max))
  
  # Si la température est en Kelvin
  data <- data |> 
    dplyr::mutate(Temperature_C = dplyr::if_else(Temperature_unit == "K", Temperature - 273.15, Temperature),
           Temperature_min_C = dplyr::if_else(Temperature_unit == "K", Temperature_min - 273.15, Temperature_min),
           Temperature_max_C = dplyr::if_else(Temperature_unit == "K", Temperature_max - 273.15, Temperature_max),
           Temperature_K = dplyr::if_else(Temperature_unit == "K", Temperature, Temperature + 273.15),
           Temperature_min_K = dplyr::if_else(Temperature_unit == "K", Temperature_min, Temperature_min + 273.15),
           Temperature_max_K = dplyr:: if_else(Temperature_unit == "K", Temperature_max, Temperature_max + 273.15))
  data <- data %>%
    dplyr::mutate(T_algo_K = dplyr::if_else(!is.na(Temperature_K), Temperature_K, (Temperature_max_K + Temperature_min_K) / 2)) %>% 
    dplyr::select(-c(Temperature_K, Temperature_min_K, Temperature_max_K))
  
  data <- data %>%
    dplyr::mutate(PAR_algo = dplyr::if_else(!is.na(PAR), PAR, (PAR_max + PAR_min) / 2))
  
  return(data)
}



select_std_or_standardisable <- function(data, sp_storing) {
  
  data <- data |> 
    dplyr::left_join(sp_storing, by = c('gragg' = 'gragg' )) |>  dplyr::mutate(Taxon = stringr::str_replace(Taxon, " ", "_"))
  ## verifier que sp storing a toutes les infos
  
  # Identifier les spcode.agg qui ne sont pas dans sp_storing
  missing_spcodes <- data |>
    dplyr::anti_join(sp_storing, by = c('gragg' = 'gragg')) |>
    dplyr::distinct(gragg)
  
  # Vérifier s'il y a des spcode.agg manquants et afficher un message
  if (nrow(missing_spcodes) > 0) {
    cat("Storing species file must be updated. The following spcode.agg here gragg are missing:\n")
    print(missing_spcodes)
  }
    
    ##
  data_std <- data |> dplyr::filter(Standardized == "true")
  
  data<- data |> dplyr::filter(Standardized == "false") |> dplyr::filter(Standardization_algo_ref != "T80" & Standardization_algo_ref != "T91" & Standardization_algo_ref != "S97") 
  
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



select_std_or_standardisable_2 <- function(data) {
  
  data_std <- data |> dplyr::filter(Standardized == "true")|> dplyr::filter(Standardization_algo_ref != "T80" & Standardization_algo_ref != "T91" & Standardization_algo_ref != "S97") 
  
  data<- data |> dplyr::filter(Standardized == "false") 
  
  data <- data |> dplyr::filter(
      
 Compound == "monoterpenes" &  ((Temperature != "NA" & PAR != "NA") | 
                                                          
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
  
  print(unique(data$Emission_unit_comp))
  return(data)
}




create_venn_diagram <- function(data) {
  # Subset the data for monoterpenes and isoprene
  monoterpenes <- dplyr::filter(data, Compound %in% c("monoterpenes"))
  isoprene <- dplyr::filter(data, Compound %in% c("isoprene"))
  
  # Remove NA values from spcode
  monoterpenes <- na.omit(monoterpenes$gragg)
  isoprene <- na.omit(isoprene$gragg)
  
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
    "Castelporziano" = c("Castelporziano NA", "Castelporziano Castelporziano", "Castelporziano Santo Querico", "Castelporziano Santo Quercio Site"),
    "Freiburg" = c("Freiburg Hartheimer Wald", "NA Hartheimer, Near Freiburg"),
    "Madrid" = c("Madrid Tres Cantos", "Madrid NA"),
    "Julich" = c("Julich Julich Research Center", "Julich NA"),
    "Hyytiala" = c("NA Hyytiala", "Hyytiala Station For Measuring Forest Ecosystem–Atmosphere Relations", "Hyytiala NA"),
    "Arad" = c("Arad Lipova Forest", "Arad NA"), "Donon" = c("Donon Vosges Forest" , "Donon Vosges"), 
    "Observatoire de haute Provence, O3HP" = c(" O3hp, Observatoire De Haute Provence" , "Marseille O3hp, Observatoire De Haute Provence"), "Terrassa" = c("Terrassa " , "Barcelona Terrassa")
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


rename_coumpounds <- function(data){
  
  
  data <- data |> 
    dplyr::mutate(Compound = stringr::str_replace(Compound, "oxygenated-monoterpenes", "ox-monoterpenes"),
                  Compound = stringr::str_replace(Compound, "oxygenated monoterpenes", "ox-monoterpenes"), 
                  Compound = stringr::str_replace(Compound, "oxygenated sesquiterpenes", "ox-sesquiterpenes")
                  )
  
  return(data)
}
# Définition de la fonction
std_iso_G93 <- function(L, T, E) {
  # Définition des constantes
  alpha <- 0.0027
  CL1 <- 1.066
  CT1 <- 95000
  CT2 <- 230000
  TM <- 314
  R <- 8.314
  TS <- 303
  
  # Calcul des termes CL et CT
  CL <- CL1 * alpha * L / sqrt(1 + alpha^2 * L^2)
  CT <- exp(CT1 * (T - TS) / (R * TS * T)) / (1 + exp(CT2 * (T - TM) / (R * TS * T)))
  
  # Calcul de ES
  ES <- E / (CL * CT)
  
  return(ES)
}

# Définition de la fonction
std_mono_G93 <- function(T, E, beta = 0.09) {
  # Définition des constantes
  TS <- 303
  
  # Calcul de CT
  CT <- exp(beta * (T - TS))
  
  # Calcul de ES
  ES <- E / CT
  
  return(ES)
}

apply_standardization <- function(data) {
  # Appliquer l'algorithme pour l'isoprène si Standardized est différent de "true"
  data <- data |> 
    dplyr::mutate(ES_iso_G93 = dplyr::if_else(Standardized != "true", std_iso_G93(L = PAR_algo, T = T_algo_K, E = Emission), NA_real_))
  
  # Appliquer l'algorithme pour les monoterpènes si Standardized est différent de "true"
  data <- data |> 
    dplyr::mutate(ES_mono_G93 = dplyr::if_else(Standardized != "true", std_mono_G93(T = T_algo_K, E = Emission), NA_real_))
  
  return(data)
}



standardisation <- function(data){
  data <- apply_standardization(data)
  
  # 1- Keep Emission when it was already standardized
  # 2- Take ES_isi_93 by default, if coumpound is monoterpenes and storing species take mono_93
  
  data |>  dplyr:: mutate(EF = dplyr::case_when(
    Standardized == "true" & !is.na(Emission) ~ Emission,
    #Stockage == "oui" & Compound =="monoterpenes" & !is.na(ES_mono_G93) ~ ES_mono_G93,
    TRUE ~ ES_iso_G93
  )) 
  
}


count_available <- function (data, minimum_nb_origin_pop){
  #espèces =Taxon pour les quelles 3 records = 3 distinct Origin_pop pour compound =  isoprene et 3 records pour monoterpenes

  species_counts <- data |> 
    dplyr::group_by(Taxon, gragg, Compound) |> 
    dplyr::summarise(distinct_origins = dplyr::n_distinct(Origin_pop), nb_trees = sum(as.numeric(Individual_replicates), na.omit = TRUE), nb_entries = dplyr::n(), .groups = 'drop')
  
  # Filter species with at least 'minimum_nb_origin_pop' distinct Origin_pop for each compound
  available_species <- species_counts |>
    dplyr::group_by(Taxon, gragg) |>
    dplyr::filter(all(distinct_origins >= minimum_nb_origin_pop)) |>
    dplyr::ungroup()
  
  # Pivot the data to have separate columns for isoprene and monoterpenes
  available_species <- available_species |>
    tidyr::pivot_wider(names_from = Compound, values_from = c(distinct_origins, nb_trees, nb_entries) , values_fill = list(distinct_origins = 0))
  
  
  return(available_species)
  
}

