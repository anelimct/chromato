sum_df_30_25draw <- readRDS("C:/Users/mchateau/Nextcloud2/THESE/chromato/data/article1/sum_df_30_25draw.rds")
sum_df_35_25draw <- readRDS("C:/Users/mchateau/Nextcloud2/THESE/chromato/data/article1/sum_df_35_25draw.rds")
sum_df_40_25draw <- readRDS("C:/Users/mchateau/Nextcloud2/THESE/chromato/data/article1/sum_df_40_25draw.rds")



library(dplyr)

metrics <- c("Richness", "functional_dispersion", "originality",
             "regularity", "dimensionality")

moyenne_par_fenetre <- function(df) {
  df |>
    group_by(window_id, median_log_sum) |>
    summarise(across(all_of(metrics), \(x) mean(x, na.rm = TRUE)),
              .groups = "drop") |>
    arrange(median_log_sum)
}

sum_df_30_25draw_mean <- moyenne_par_fenetre(sum_df_30_25draw)
sum_df_35_25draw_mean <- moyenne_par_fenetre(sum_df_35_25draw)
sum_df_40_25draw_mean <- moyenne_par_fenetre(sum_df_40_25draw)





library(dplyr)

metrics <- c("Richness", "functional_dispersion", "originality",
             "regularity", "dimensionality")

mediane_par_fenetre <- function(df) {
  df |>
    group_by(window_id, median_log_sum) |>
    summarise(across(all_of(metrics), \(x) median(x, na.rm = TRUE)),
              .groups = "drop") |>
    arrange(median_log_sum)
}

sum_df_30_25draw_median <- mediane_par_fenetre(sum_df_30_25draw)
sum_df_35_25draw_median <- mediane_par_fenetre(sum_df_35_25draw)
sum_df_40_25draw_median <- mediane_par_fenetre(sum_df_40_25draw)



library(dplyr)

metrics <- c("Richness", "functional_dispersion", "originality",
             "regularity", "dimensionality")

q05_par_fenetre <- function(df) {
  df |>
    group_by(window_id, median_log_sum) |>
    summarise(across(all_of(metrics), \(x) quantile(x, probs = 0.05, na.rm = TRUE, names = FALSE)),
              .groups = "drop") |>
    arrange(median_log_sum)
}

sum_df_30_25draw_q05 <- q05_par_fenetre(sum_df_30_25draw)
sum_df_35_25draw_q05 <- q05_par_fenetre(sum_df_35_25draw)
sum_df_40_25draw_q05 <- q05_par_fenetre(sum_df_40_25draw)

library(dplyr)

metrics <- c("Richness", "functional_dispersion", "originality",
             "regularity", "dimensionality")

q95_par_fenetre <- function(df) {
  df |>
    group_by(window_id, median_log_sum) |>
    summarise(across(all_of(metrics), \(x) quantile(x, probs = 0.95, na.rm = TRUE, names = FALSE)),
              .groups = "drop") |>
    arrange(median_log_sum)
}

sum_df_30_25draw_q95 <- q95_par_fenetre(sum_df_30_25draw)
sum_df_35_25draw_q95 <- q95_par_fenetre(sum_df_35_25draw)
sum_df_40_25draw_q95 <- q95_par_fenetre(sum_df_40_25draw)

