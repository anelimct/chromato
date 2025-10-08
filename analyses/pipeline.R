# Created by use_targets().
# Follow the comments below to fill in this target script.
# Then follow the manual to check and run the pipeline:
#   https://books.ropensci.org/targets/walkthrough.html#inspect-the-pipeline

# Load packages required to define the pipeline:
library(targets)
library(tarchetypes)
library(ggplot2)
library(lubridate)
# library(tarchetypes) # Load other packages as needed.

# Set target options:
tar_option_set(
  packages = c("tibble") # Packages that your targets need for their tasks.
  # format = "qs", # Optionally set the default storage format. qs is fast.
  #
  # Pipelines that take a long time to run may benefit from
  # optional distributed computing. To use this capability
  # in tar_make(), supply a {crew} controller
  # as discussed at https://books.ropensci.org/targets/crew.html.
  # Choose a controller that suits your needs. For example, the following
  # sets a controller that scales up to a maximum of two workers
  # which run as local R processes. Each worker launches when there is work
  # to do and exits if 60 seconds pass with no tasks to run.
  #
  #   controller = crew::crew_controller_local(workers = 2, seconds_idle = 60)
  #
  # Alternatively, if you want workers to run on a high-performance computing
  # cluster, select a controller from the {crew.cluster} package.
  # For the cloud, see plugin packages like {crew.aws.batch}.
  # The following example is a controller for Sun Grid Engine (SGE).
  # 
  #   controller = crew.cluster::crew_controller_sge(
  #     # Number of workers that the pipeline can scale up to:
  #     workers = 10,
  #     # It is recommended to set an idle time so workers can shut themselves
  #     # down if they are not running tasks.
  #     seconds_idle = 120,
  #     # Many clusters install R as an environment module, and you can load it
  #     # with the script_lines argument. To select a specific verison of R,
  #     # you may need to include a version string, e.g. "module load R/4.3.2".
  #     # Check with your system administrator if you are unsure.
  #     script_lines = "module load R"
  #   )
  #
  # Set other options as needed.
)

# Run the R scripts in the R/ folder with your custom functions:
tar_source()
# tar_source("other_functions.R") # Source other scripts as needed.

# Replace the target list below with your own:
list(
  
  tar_target(alcanes_samples_file, here::here("data", "alcanes_samples.csv" ), format = "file"),
  tar_target(alcanes_samples,utils::read.csv(alcanes_samples_file, sep = ";")),
  tar_target(renamed_alcanes, rename_GC_files(alcanes_samples, "alcanes") |> dplyr::filter(exploitables == "oui") |> dplyr::mutate(dplyr::across(c("octane", "nonane", "decane", "undecane", "dodecane", 
     
                                                                                                                                                                                                                                                                                          "tridecane", "tetradecane", "pentadecane", "hexadecane"), as.numeric)) |> imput_missing_alcanes() ) ,
  tar_target(calib_quanti_file, here::here("data", "calib_quanti.csv" ), format = "file"),
  tar_target(calib_quanti,utils::read.csv(calib_quanti_file, sep = ";")),
  
  
  tar_target(bvocs_samples_file, here::here("data", "BVOCs_samples.csv" ), format = "file"),
  tar_target(bvocs_samples_raw,utils::read.csv(bvocs_samples_file, sep = ";")),
  
  
  ##Ranger les chromato par batch 
  tar_target(create_files_batch_2023, organize_gc_files_by_batch(subset_2023, "2023" )),
  tar_target(create_files_batch_2024, organize_gc_files_by_batch(subset_2024, "2024" )),
  tar_target(create_files_batch_2025, organize_gc_files_by_batch(subset_2025, "2025" )),
  
  
  tar_target(woodiv_species_file, here::here("data", "WOODIV_DB_release_v2", "SPECIES",  "WOODIV_v2_Species_code.csv"), format = "file"), 
  tar_target(woodiv_species, utils::read.csv(woodiv_species_file)),
  
<<<<<<< HEAD
  tar_target(terrain_2025_file, here::here("data", "terrain_2025.csv"), format = "file"), 
  tar_target(terrain_2025, utils::read.csv(terrain_2025_file, sep =";")),
  
  tar_target(tree, ape::read.tree( paste0( here::here("data", "WOODIV") , "/INTEGRADIV_phylogeny_trees.tre" ))), 
=======
  tar_target(tree, ape::read.tree( paste0( here::here("data", "WOODIV_DB_release_v2", "PHYLOGENY") , "/WOODIV_v2_Phylogeny_gragg.tree" ))), 
>>>>>>> temp-branch
  
  #liste des paradise reports
   #mono
  tar_target( paradise_reports_files, list.files(here::here("data", "paradise_reports"), pattern = "\\.xlsx$")),
  tar_target( paradise_reports_list , read_paradise_4( paradise_reports_files, library_CAS_RI)),
  
  #iso
  tar_target( paradise_reports_iso_files, list.files(here::here("data", "paradise_reports", "iso"), pattern = "\\.xlsx$")),
  tar_target( paradise_reports_iso_list , read_paradise_4( paradise_reports_iso_files, library_CAS_RI, option = "iso/")),
  
  #liste des paradise reports calibrations iso et mono, pour avoir à la fois tous les noms des fichiers CDF calib mais aussi la courbe de calibration a partir d'un batch sans samples = TRUE calibration curve to compare with batches calibration curves
  tar_target( paradise_reports_calib_files, list.files(here::here("data", "paradise_reports", "calib"), pattern = "\\.xlsx$")),
  tar_target( paradise_reports_calib_list , read_paradise_4( paradise_reports_calib_files, library_CAS_RI, option = "calib/")),
  
  ## Récupérer les données des iButtuns
  
  tar_target(ibuttons_files, list.files ( here::here("data", "iButton_files"), full.names = TRUE), format = "file"), 
  tar_target(ibutton_table_T, make_ibuttons_table(ibuttons_files, "T")), 
  tar_target(ibutton_table_H, make_ibuttons_table(ibuttons_files, "H")),
  tar_target(bvocs_ibutton_values_T, export_ibuttons_data(bvocs_samples_raw, ibutton_table_T, "T")),
  tar_target(bvocs_ibutton_values_H, export_ibuttons_data(bvocs_samples_raw, ibutton_table_H, "H")), 
  tar_target(bvocs_samples_ibuttons_values, dplyr::left_join(bvocs_ibutton_values_T, bvocs_ibutton_values_H) |>  create_new_par_column() ) ,

             
  ## Trier les chromato, subset is done on desorption date ? or sample date ? 
  
  tar_target(subset_2023, subset_year(bvocs_samples_ibuttons_values, "2023", renamed_alcanes) |> save_plot_ibuttons ("2023", plot_in_out) |> save_plot_corr_T ("2023", plot_corr_Tin_Tout)),
  tar_target(subset_2024, subset_year(bvocs_samples_ibuttons_values, "2024", renamed_alcanes) |>  save_plot_ibuttons( "2024", plot_in_out) |> save_plot_corr_T ( "2024", plot_corr_Tin_Tout) ), 
  tar_target(subset_2025, subset_year(bvocs_samples_ibuttons_values, "2025", renamed_alcanes) |>  save_plot_ibuttons( "2025", plot_in_out) |> save_plot_corr_T ( "2025", plot_corr_Tin_Tout) ),
  
  tar_target(check_chromato_2024, chromato_file_check(subset_2024, "2024")), 
  tar_target(check_chromato_2023, chromato_file_check(subset_2023, "2023")),
  tar_target(check_chromato_2025, chromato_file_check(subset_2025, "2025")), 
  
  tarchetypes::tar_quarto(report, "01_presentation_batches.qmd"),
  #tarchetypes::tar_quarto(report, "02_presentation_field.qmd"),
  
  ##All bvocs samples lines with values i buttons as a list, closest alcane, batch, if we can read or not the sample (col iso, col mono TRUE or FALSE)
  tar_target(bvocs_samples, rbind(subset_2023, subset_2024, subset_2025) |> dplyr::select(- Time_T_in, -Time_T_out, -Time_H_in, -Time_H_out, -values_H_out, -start_prelevement, -end_prelevement ) |> var_paradise( paradise_reports_list,paradise_reports_iso_list, calib_quanti)|> paired_samples()),
  
  ##Filter according to a set of conditions
  tar_target(failed_samples, filter_out_samples(bvocs_samples, 43, 42, type = "failed", "mono")),#43 Tmax, mean 42
  tar_target(valid_samples_mono,  filter_out_samples(bvocs_samples, 43, 42, type = "keep", "mono") |> species_aggregation( woodiv_species, "field")),
  tar_target(valid_samples_iso,  filter_out_samples(bvocs_samples, 43, 42, type = "keep", "iso") |> species_aggregation( woodiv_species, "field")),
  
  tar_target(summary_field, summarize_field (valid_samples_iso, valid_samples_mono)),
  
  
  ##Lire les articles BD_litt
  
  tar_target(articles, list.files(here::here("data", "TREEVOCS_data", "TREEVOCS_data_shiny"), full.names = TRUE), format = "file"), 
  
  tar_target(DB_bvocs, {
    files <- list()
    for(file in articles) {
      df <- read_excel_articles(file)
      files[[file]] <- df
    }
    
    # Check consistency
    col_counts <- sapply(files, ncol)
    expected_cols <- names(files[[1]])  # Use first file as reference
    
    problematic <- c()
    for(i in seq_along(files)) {
      if(!identical(names(files[[i]]), expected_cols)) {
        problematic <- c(problematic, articles[i])
        cat("PROBLEM with:", articles[i], "\n")
        cat(" Issue in columns:  ", paste(setdiff(names(files[[i]]), expected_cols ), collapse=", "), "\n")
        cat("  Column count:", ncol(files[[i]]), "\n\n")
      }
    }
    
    if(length(problematic) > 0) {
      stop("Column mismatch in files: ", paste(problematic, collapse=", "))
    }
    
    do.call(rbind, files)
  } |> species_aggregation(woodiv_species, "litt") ), 

  
  tar_target(DB_bvocs_filtered, DB_bvocs |>  select_iso_mono() |> numeric_emissions_g_h() |> dplyr::filter(Emission_unit_leaf == "g" | Emission == 0) |> select_std_or_standardisable_2()  |>  select_months( "05", "07")  |>  select_temp_and_par(42, 20, 1500, 500) |> 
               dplyr::mutate(
                 Country = clean_country(Country),
                 Origin_city = clean_city_locality(Origin_city),
                 Origin_locality = clean_city_locality(Origin_locality),
                 Origin_pop = paste(Origin_city, Origin_locality, sep = " "),
                 Origin_pop = ifelse(Origin_pop == "NA NA", Ref_ID_WoS, Origin_pop)
               ) |> create_population_variable ()),
  #tar_target(DB_bvocs_ES, standardisation (DB_bvocs_filtered) |>  boxplot_EF(tree, field_EF)),
  
  #tar_target(summary_DB, count_available (DB_bvocs_ES, 1)),
  
  ## Library CAS

  tar_target(RI_file, here::here("data", "library_cas_ri.xlsx" ), format = "file"),
  tar_target(RI_th,  readxl::read_xlsx(RI_file, sheet = 1)), 
  tar_files(CAS_files, list.files(here::here("data", "web_requests_CAS"), full.names = TRUE)
  ), 
  
  tar_target(library_CAS, create_library(CAS_files), pattern = map(CAS_files)),
  tar_target(library_CAS_RI, update_lib(library_CAS, RI_th)),
  


#   ## working on paradise_files
#
#
  tar_target(RI_exp ,compute_retention_index( renamed_alcanes, bvocs_samples, paradise_reports_list, calib_quanti )),
  tar_target(table_calib_mono_btw_session, compare_calib_btw_reports(calib_quanti, paradise_reports_list, paradise_reports_calib_list, "mono") |>   plot_calib_btw_session( calib_quanti, "mono", library_CAS_RI)),
<<<<<<< HEAD
  
  tar_target(table_calib_iso_btw_session, compare_calib_btw_reports(calib_quanti, paradise_reports_iso_list, paradise_reports_calib_list, "iso")|>   plot_calib_btw_session( calib_quanti, "iso", library_CAS_RI)),
  
  tar_target( table_blanks_list, create_list_dataframes_blanks(bvocs_samples) |> sort_list_blanks_tables()),
  
  tar_target(paradise_grouped_blanks_iso, paradise_reports_grouped_blanks(paradise_reports_iso_list, table_blanks_list, iso = TRUE)),
  
  tar_target(paradise_reports_sbtr_blanks_iso_list, subtract_blanks_from_samples(paradise_reports_iso_list, paradise_grouped_blanks_iso, calib_quanti)), 
  
  tar_target(paradise_reports_iso_quanti_list, area_to_quanti(paradise_reports_sbtr_blanks_iso_list, calib_quanti, table_calib_iso_btw_session)), 
  
  tar_target(paradise_reports_iso_ER_list, compute_ER (paradise_reports_iso_quanti_list, bvocs_samples, calib_quanti)), 
  
  tarchetypes::tar_quarto(report, "01_presentation_batches.qmd"),
  
  
  ##SPATIAL maps
  tar_target(working_file, {
    utils::read.csv(paste0(here::here("data", "WOODIV_DB_release_v2", "OCCURRENCE"), "/WOODIV_v2_Occurrence_data.csv")) |> dplyr::left_join(woodiv_species, by = "spcode")
  }), 
  tar_target(WOODIV_grid, {
    WOODIV_grid <- sf::st_read(paste0(here::here("data", "WOODIV_DB_release_v2", "SPATIAL", "WOODIV_v2_Grid_epsg3035"), "/WOODIV_v2_Grid_epsg3035.shp"))
    WOODIV_grid <- WOODIV_grid[WOODIV_grid$idgrid %in% unique(working_file$idgrid), ]
  }),
  tar_target(WOODIV_shape, {
    sf::st_read(paste0(here::here("data", "WOODIV_DB_release_v2", "SPATIAL", "WOODIV_v2_Shape_epsg3035"), "/WOODIV_v2_Shape_epsg3035.shp"))
  }),
 
  
  tar_target(summary_all , ranking_species(working_file) |>  dplyr::left_join(summary_DB, by = c('gragg' = 'gragg')) |>  dplyr::left_join(summary_field,  by = c('gragg' = 'gragg')) |> 
               process_summary_data(working_file)  |>  plot_hist_ranking("all_distinct_origins_isoprene", 20, "Effort échantilonnage isoprène pour les espèces les plus communes ") |> plot_hist_ranking("all_distinct_origins_monoterpenes", 20, "Effort échantilonnage monoterpènes pour les espèces les plus communes ") |>  plot_tree_effort_ech(tree) |>  plot_hist_ranking_cumul (terrain_2025) 
             |> tidy_summary_all()),
              
  tar_target(completeness , compute_completeness(WOODIV_grid, working_file, summary_all, 1) |> map_et_plot_completness(WOODIV_shape))
  
  #How important are species to sample to have max completness
=======
  tar_target(table_calib_iso_btw_session, compare_calib_btw_reports(calib_quanti, paradise_reports_iso_list, paradise_reports_calib_list, "iso") |>   plot_calib_btw_session( calib_quanti, "iso", library_CAS_RI)),
>>>>>>> temp-branch

#   ##ISOPRENE
#
  tar_target(paradise_reports_sbtr_blanks_iso_list, subtract_blanks_from_samples(paradise_reports_iso_list, bvocs_samples, calib_quanti)),
  tar_target(paradise_reports_iso_quanti_list, area_to_quanti_iso(paradise_reports_sbtr_blanks_iso_list, calib_quanti, table_calib_iso_btw_session)),
  tar_target(paradise_reports_iso_ER, compute_ER (paradise_reports_iso_quanti_list, bvocs_samples, calib_quanti) |>  mark_values(paradise_reports_iso_list, lod_3x = 193500, lod_10x = 645000) |>  sum_isoprene_across_reports()),
#
#
#   ##MONO/SESQUI
  tar_target(paradise_reports_sbtr_blanks_mono_list, subtract_blanks_from_samples(paradise_reports_list, bvocs_samples, calib_quanti)),
  tar_target(paradise_reports_mono_quanti_list, area_to_quanti_mono(paradise_reports_sbtr_blanks_mono_list, calib_quanti, table_calib_mono_btw_session)),
  tar_target(paradise_reports_mono_ER, compute_ER (paradise_reports_mono_quanti_list, bvocs_samples, calib_quanti) |>  mark_values(paradise_reports_list, lod_3x = 193500, lod_10x = 645000)|>  lapply(chemodiv::NPCTable) |>  sum_terpenoids_across_reports()),
  tar_target(field_EF ,  merge_datasets (paradise_reports_mono_ER, paradise_reports_iso_ER, valid_samples_iso) |>  species_aggregation(woodiv_species, "field_"))
#
#
#
#
#   ##SPATIAL maps
#   tar_target(working_file, {
#     utils::read.csv(paste0(here::here("data", "WOODIV_DB_release_v2", "OCCURRENCE"), "/WOODIV_v2_Occurrence_data.csv")) |> dplyr::left_join(woodiv_species, by = "spcode")
#   }),
#   tar_target(WOODIV_grid, {
#     WOODIV_grid <- sf::st_read(paste0(here::here("data", "WOODIV_DB_release_v2", "SPATIAL", "WOODIV_v2_Grid_epsg3035"), "/WOODIV_v2_Grid_epsg3035.shp"))
#     WOODIV_grid <- WOODIV_grid[WOODIV_grid$idgrid %in% unique(working_file$idgrid), ]
#   }),
#   tar_target(WOODIV_shape, {
#     sf::st_read(paste0(here::here("data", "WOODIV_DB_release_v2", "SPATIAL", "WOODIV_v2_Shape_epsg3035"), "/WOODIV_v2_Shape_epsg3035.shp"))
#   }),
#
#
#   tar_target(summary_all , ranking_species(working_file) |>  dplyr::left_join(summary_DB, by = c('gragg' = 'gragg')) |>  dplyr::left_join(summary_field,  by = c('gragg' = 'gragg')) |>
#                process_summary_data(working_file)  |>  plot_hist_ranking("all_distinct_origins_isoprene", 20, "Effort échantilonnage isoprène pour les espèces les plus communes ") |> plot_hist_ranking("all_distinct_origins_monoterpenes", 20, "Effort échantilonnage monoterpènes pour les espèces les plus communes ") |>  plot_tree_effort_ech(tree) |>  plot_hist_ranking_cumul (terrain_2025)
#              |> tidy_summary_all(woodiv_species)),
#
#   tar_target(completeness , compute_completeness(WOODIV_grid, working_file, summary_all, 1) |> map_et_plot_completness(WOODIV_shape))
#
#   #How important are species to sample to have max completness
#
  )
