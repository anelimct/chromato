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
  
  #tar_target(calib_samples_file, here::here("data", "calib_samples.csv" ), format = "file"),
  #tar_target(calib_samples,utils::read.csv(calib_samples_file, sep = ";")),
  #tar_target(renamed_calib, rename_GC_files(calib_samples, "calib")),
  
  tar_target(calib_quanti_file, here::here("data", "calib_quanti.csv" ), format = "file"),
  tar_target(calib_quanti,utils::read.csv(calib_quanti_file, sep = ";")),
  
  
  tar_target(bvocs_samples_file, here::here("data", "BVOCs_samples.csv" ), format = "file"),
  tar_target(bvocs_samples_raw,utils::read.csv(bvocs_samples_file, sep = ";")),
  
  tar_target(storing_species_file, here::here("data", "splist_storing_AS.csv"), format = "file"), 
  tar_target(storing_species, utils::read.csv(storing_species_file, sep = ";")),
  
  tar_target(woodiv_species_file, here::here("data", "WOODIV", "WOODIV_Species.csv")), 
  tar_target(woodiv_species, utils::read.csv(woodiv_species_file)),
  
  #liste des paradise reports
  
  tar_target( paradise_reports_files, list.files(here::here("data", "paradise_reports"), pattern = "\\.xlsx$")),
  tar_target( paradise_reports_list , read_paradise_4( paradise_reports_files, library_CAS_RI)),
  
  #liste des paradise reports calibrations
  
  tar_target( paradise_reports_calib_files, list.files(here::here("data", "paradise_reports", "calib"), pattern = "\\.xlsx$")),
  tar_target( paradise_reports_calib_list , read_paradise_4( paradise_reports_calib_files, library_CAS_RI, option = "calib/")),
  
  ## Récupérer les données des iButtuns
  
  tar_target(ibuttons_files, list.files ( here::here("data", "iButton_files"), full.names = TRUE), format = "file"), 
  tar_target(ibutton_table_T, make_ibuttons_table(ibuttons_files, "T")), 
  tar_target(ibutton_table_H, make_ibuttons_table(ibuttons_files, "H")),
  tar_target(bvocs_ibutton_values_T, export_ibuttons_data(bvocs_samples_raw, ibutton_table_T, "T")),
  tar_target(bvocs_ibutton_values_H, export_ibuttons_data(bvocs_samples_raw, ibutton_table_H, "H")), 
  tar_target(bvocs_samples_ibuttons_values, dplyr::left_join(bvocs_ibutton_values_T, bvocs_ibutton_values_H)),
  ##Lire les articles BD_litt
  
  tar_target(articles, list.files(here::here("data", "TREEVOCS_data", "TREEVOCS_data_extraction_R_edit_v_stage"), full.names = TRUE), format = "file"), 
  tar_target(DB_bvocs, {
    files <- lapply(articles, read_excel_articles)
    do.call(rbind, files)
  } |> species_aggregation(woodiv_species)|> select_iso_mono() |> numeric_emissions_g_h() |> dplyr::filter(Emission_unit_leaf == "g" | Emission == 0) ),
  
  tar_target(DB_bvocs_filtered, select_std_or_standardisable(DB_bvocs, storing_species)  |>  select_months( "05", "07")  |>  select_temp_and_par(42, 20, 1500, 500) |> 
               dplyr::mutate(
                 Country = clean_country(Country),
                 Origin_city = clean_city_locality(Origin_city),
                 Origin_locality = clean_city_locality(Origin_locality),
                 Origin_pop = paste(Origin_city, Origin_locality, sep = " "),
                 Origin_pop = ifelse(Origin_pop == "NA NA", Ref_ID_WoS, Origin_pop)
               ) |> create_population_variable ()),
             
  ## Trier les chromato, subset is done on desorption date
  
  tar_target(subset_2023, subset_year(bvocs_samples_ibuttons_values, "2023", renamed_alcanes)),
  tar_target(subset_2024, subset_year(bvocs_samples_ibuttons_values, "2024", renamed_alcanes)), 
  
  tar_target(check_chromato_2024, chromato_file_check(subset_2024, "2024")), 
  tar_target(check_chromato_2023, chromato_file_check(subset_2023, "2023")),
  
  tar_target(chronologie_2023, chronologie(subset_2023)), 
  tar_target(chronologie_2024, chronologie(subset_2024)), 
  
  tar_target(bvocs_samples, rbind(subset_2023, subset_2024) |> dplyr::select(- Time_T_in, -Time_T_out, -Time_H_in, -Time_H_out, -values_H_out, -start_prelevement, -end_prelevement ) |> var_paradise( paradise_reports_list, calib_quanti ) |> paired_samples()),
  

  tar_target(failed_samples, filter_out_samples(bvocs_samples, 42, 42)), #42 Tmax mean 40
  
  
  
  ##Ranger les chromato par batch avec leurs alcanes correspondants
  tar_target(create_files_batch_2023, organize_gc_files_by_batch(subset_2023, "2023" )),
  tar_target(create_files_batch_2024, organize_gc_files_by_batch(subset_2024, "2024" )),
  
  tar_target(plot_ibuttons_2024, save_plot_ibuttons (subset_2024, "2024", plot_in_out)),
  tar_target(plot_corr_T_2024, save_plot_corr_T (subset_2024, "2024", plot_corr_Tin_Tout)),
  tar_target(plot_ibuttons_2023, save_plot_ibuttons (subset_2023, "2023", plot_in_out)), 
  tar_target(plot_corr_T_2023, save_plot_corr_T (subset_2023, "2023", plot_corr_Tin_Tout)),
 
  ## Library CAS

  tar_target(RI_file, here::here("data", "library_cas_ri.xlsx" ), format = "file"),
  tar_target(RI_th,  readxl::read_xlsx(RI_file, sheet = 1)), 
  tar_files(
    CAS_files,
    list.files(here::here("data", "web_requests_CAS"), full.names = TRUE)
  ), 
  
  tar_target(library_CAS, create_library(CAS_files), pattern = map(CAS_files)),
  tar_target(library_CAS_RI, update_lib(library_CAS, RI_th)),
  

  
  ## working on paradise_files
  
  
  tar_target(RI_exp ,compute_retention_index( renamed_alcanes, bvocs_samples, paradise_reports_list, calib_quanti )),
  
  tarchetypes::tar_quarto(report, "01_presentation_batches.qmd")

)
