plot_compare_standardisation <- function(compounds_table_standardized, valid_samples_mono, times_compound_sp, bvocs_samples) {
  #mettre aussi bvocs samples pour avoir les T est les PAR pour voir si ça a un impact sur la standardisation
  #pour les all samples il y a aussi all compounds meme les uniques pour cette espèces
  # spagg_to_keep <- unique(times_compound_sp$spagg)
  # pour les all_samples il y ameme toutes les espèces que l'on a fait pas que celles pour l'article de screening,
  #comme ça ça nous permet aussi de tester ce qui sera dans la DB total et ce dont on se servira pour la suite 
  
  all_samples_T <- compounds_table_standardized[[3]]
  all_samples_T_L <- compounds_table_standardized[[2]]
  all_samples_T_bourtsou <- compounds_table_standardized[[4]]
  
  cols_PAR <- c("PAR.début.1", "PAR.début.2", "PAR.milieu.1",     
                "PAR.milieu.2", "PAR.milieu.3", "PAR.milieu.4", "PAR.milieu.5", "PAR.milieu.6", "PAR.fin.1", "PAR.fin.2")
  
  bvocs_samples<- bvocs_samples |> dplyr::mutate(across(all_of(cols_PAR), as.numeric)) |>  
    dplyr::mutate(PAR_algo = rowMeans(dplyr::across(all_of(cols_PAR)), na.rm = TRUE)) |> 
    dplyr::mutate(mean_T = purrr::map_dbl(values_T_in, ~ mean(.x))) |>  dplyr::filter(Taxon != "BLANC") |>  dplyr::mutate(sample = stringr::str_c(stringr::str_to_lower(ID), ".cdf")) |>  dplyr::select(sample, PAR_algo, mean_T)
    ## la il faudrait avant faire la colonne mean PAR et mean T qui je en sais pas pourquoi je la refait a chaque fois et je en l'enregistre jamais dans le target c'est insuportable

  all_samples_T_sum <- sum_byclass(compounds_table_standardized[[3]]) |>  tidyr::pivot_longer(cols = 2:305, names_to = "sample", values_to = "Emission") |> dplyr::left_join(bvocs_samples)
  all_samples_T_L_sum <- sum_byclass(compounds_table_standardized[[2]])|>  tidyr::pivot_longer(cols = 2:305, names_to = "sample", values_to = "Emission")|> dplyr::left_join(bvocs_samples)
  all_samples_T_bourtsou_sum <- sum_byclass(compounds_table_standardized[[4]])|>  tidyr::pivot_longer(cols = 2:305, names_to = "sample", values_to = "Emission")|> dplyr::left_join(bvocs_samples)
  
  # Créer le plot horizontal
  plot_1 <- ggplot(all_samples_T_sum, aes(x = PAR_algo, y = Emission, color = class)) +
    geom_point(position = "jitter", alpha = 0.8) +
    scale_fill_manual(
      values = c("Monoterpenes" = "steelblue", 
                 "Oxygenated-monoterpenes" = "darkorange")
    ) +
    labs(
      title = "Emission standardisées en fonction du PAR de l'échantillon",
      x = "PAR",
      y = "Emission std",
      fill = "Classe"
    ) +
    theme_minimal() +
    theme(
      legend.position = "top"
    )
  
  
  
  
  plot_2 <- ggplot(all_samples_T_sum, aes(x = mean_T, y = Emission, color = class)) +
    geom_point(position = "jitter", alpha = 0.8) +
    scale_fill_manual(
      values = c("Monoterpenes" = "steelblue", 
                 "Oxygenated-monoterpenes" = "darkorange")
    ) +
    labs(
      title = "Emission standardisées en fonction de la température de l'échantillon",
      x = "Température",
      y = "Emission std ",
      fill = "Classe"
    ) +
    theme_minimal() +
    theme(
      legend.position = "top"
    )

  
  
  mean_spagg_T <- compound_mean_sp(compounds_table_standardized[[3]], valid_samples_mono, times_compound_sp, include_se = FALSE) |> dplyr::mutate(dplyr::across(5:50, as.numeric))
  mean_spagg_T_bourstou <- compound_mean_sp(compounds_table_standardized[[4]], valid_samples_mono, times_compound_sp, include_se = FALSE)|> dplyr::mutate(dplyr::across(5:50, as.numeric))
  mean_spagg_T_L <- compound_mean_sp(compounds_table_standardized[[2]], valid_samples_mono, times_compound_sp, include_se = FALSE)|> dplyr::mutate(dplyr::across(5:50, as.numeric))
  
  
  
  mean_spagg_T_sum <- sum_byclass( mean_spagg_T)
  mean_spagg_T_bourstou_sum <- sum_byclass( mean_spagg_T_bourstou)
  mean_spagg_T_L_sum <- sum_byclass( mean_spagg_T_L)
  
  
  delta_T_vs_bourtsou_sum_spagg <-  mean_spagg_T_sum
  delta_T_vs_bourtsou_sum_spagg[-1] <- mean_spagg_T_sum[-1] - mean_spagg_T_bourstou_sum[-1] 
  
  delta_T_vs_T_L_sum_spagg <-  mean_spagg_T_sum
  delta_T_vs_T_L_sum_spagg [-1] <- mean_spagg_T_sum[-1] - mean_spagg_T_L_sum[-1] 
  
  delta_T_bourtsou_sum_spagg_vs_T_L <-   mean_spagg_T_bourstou_sum
  delta_T_bourtsou_sum_spagg_vs_T_L [-1] <-  mean_spagg_T_bourstou_sum[-1] - mean_spagg_T_L_sum[-1] 
  
  
  # Préparer les données en format long
  data_long <- delta_T_vs_T_L_sum_spagg |> 
    tidyr::pivot_longer(
      cols = -class,  # Garder la colonne class comme identifiant
      names_to = "espece", 
      values_to = "mean_value"
    )
  
  # Créer le plot horizontal
  plot <- ggplot(data_long, aes(x = mean_value, y = espece, fill = class)) +
    geom_col(position = "dodge", alpha = 0.8) +
    scale_fill_manual(
      values = c("Monoterpenes" = "steelblue", 
                 "Oxygenated-monoterpenes" = "darkorange")
    ) +
    labs(
      title = "Différence entre le facteur d'emission standardisé avec G93 T versus G93 T+L",
      x = expression(Δ*"G93T - G93T+L"),
      y = "Espèces",
      fill = "Classe"
    ) +
    theme_minimal() +
    theme(
      legend.position = "top"
    )
  
  ggsave( "std_Tvs_T_L.png",   plot, path = paste0(here::here ("figures", "standardisation")), width = 16, height =  16, units = "cm", bg = "white", create.dir = TRUE )
  
  
  # Préparer les données en format long
  data_long <- delta_T_vs_bourtsou_sum_spagg |> 
    tidyr::pivot_longer(
      cols = -class,  # Garder la colonne class comme identifiant
      names_to = "espece", 
      values_to = "mean_value"
    )
  
  # Créer le plot horizontal
  plot <-ggplot(data_long, aes(x = mean_value, y = espece, fill = class)) +
    geom_col(position = "dodge", alpha = 0.8) +
    scale_fill_manual(
      values = c("Monoterpenes" = "steelblue", 
                 "Oxygenated-monoterpenes" = "darkorange")
    ) +
    labs(
      title = "Différence entre le facteur d'emission standardisé avec G93 T constente origine (0.10) vs Bourtsoukidis meta-analysis (0.14)  ",
      x = expression(Δ*"G93T - G93T bourtsoukidis"),
      y = "Espèces",
      fill = "Classe"
    ) +
    theme_minimal() +
    theme(
      legend.position = "top"
    )
  
  ggsave( "std_Tvs_T_bourtsoukis.png",   plot, path = paste0(here::here ("figures", "standardisation")), width = 16, height =  16, units = "cm", bg = "white", create.dir = TRUE )
  
  # Préparer les données en format long
  data_long <-  delta_T_bourtsou_sum_spagg_vs_T_L  |> 
    tidyr::pivot_longer(
      cols = -class,  # Garder la colonne class comme identifiant
      names_to = "espece", 
      values_to = "mean_value"
    )
  
  # Créer le plot horizontal
  plot <-ggplot(data_long, aes(x = mean_value, y = espece, fill = class)) +
    geom_col(position = "dodge", alpha = 0.8) +
    scale_fill_manual(
      values = c("Monoterpenes" = "steelblue", 
                 "Oxygenated-monoterpenes" = "darkorange")
    ) +
    labs(
      title = "Différence entre le facteur d'emission standardisé avec G93 T+L vs G93 T constante Bourtsoukidis meta-analysis (0.14)  ",
      x = expression(Δ*" G93T bourtsoukidis - G93T+L"),
      y = "Espèces",
      fill = "Classe"
    ) +
    theme_minimal() +
    theme(
      legend.position = "top"
    )
  
  ggsave( "std_T_bourtsoukis_vs_T_L.png",   plot, path = paste0(here::here ("figures", "standardisation")), width = 16, height =  16, units = "cm", bg = "white", create.dir = TRUE )
  
  
  }


sum_byclass <- function(df){
  
  numeric_cols <- 5:(ncol(df)-1)
  
  result <- df |> 
    dplyr::group_by(class) |> 
    dplyr::summarise(
      dplyr::across(
        .cols = all_of(numeric_cols),
        .fns = ~ sum(., na.rm = TRUE),
        .names = "{.col}"
      ),
      .groups = "drop"
    )
}

