
get_compound_stats <- function(compound_table){

  data_columns <- compounds_table[, 5:312] |>
    dplyr::mutate(across(everything(), as.character))

  # Recréer le dataframe avec les colonnes converties
  compounds_table_char <- cbind(compounds_table[, 1:4], data_columns)

  # Maintenant compter
  tableau_final <- compounds_table_char |>
    dplyr::rowwise() |>
    dplyr::mutate(
      tr = sum(c_across(5:312) == "tr", na.rm = TRUE),
      nd = sum(c_across(5:312) == "nd", na.rm = TRUE),
      quantifiee = sum(!is.na(c_across(5:312)) &
                         c_across(5:312) != "tr" &
                         c_across(5:312) != "nd" &
                         c_across(5:312) != "0", na.rm = TRUE)
    ) |>
    dplyr::select(compound, tr, nd, quantifiee)

}

