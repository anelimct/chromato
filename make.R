#' chromato: A Research Compendium
#' 
#' @description 
#' A paragraph providing a full description of the project and describing each 
#' step of the workflow.
#' 
#' @author Miléna Chateau \email{milena.chateau@imbe.fr}
#' 
#' @date 2024/12/11



## Install Dependencies (listed in DESCRIPTION) ----

devtools::install_deps(upgrade = "never")


## Load Project Addins (R Functions and Packages) ----

devtools::load_all(here::here())
library(targets)
library(ggplot2)


targets::tar_config_set(store = "outputs/pipeline/", script = "analyses/pipeline.R")

## Global Variables ----

## Run Project ----

# List all R scripts in a sequential order and using the following form:
# source(here::here("analyses", "script_X.R"))
