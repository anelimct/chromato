library(purrr)
library(dplyr)
library(kableExtra)

metric_names <- c(
  Richness              = "Functional richness",
  regularity            = "Regularity",
  functional_dispersion = "Functional dispersion",
  originality           = "Originality",
  dimensionality        = "Functional dimensionality"
)

get_aic <- function(aic_table, model) {
  v <- aic_table$AIC[aic_table$model == model]
  if (length(v) == 0) NA_real_ else v
}

fmt_thr <- function(thr, lo, hi, i) {
  if (length(thr) < i || is.na(thr[i])) return("--")
  sprintf("%.2f [%.2f ; %.2f]", thr[i], lo[i], hi[i])
}

tab <- imap_dfr(results, function(win, win_nm) {
  imap_dfr(win, function(res, metric) {
    best <- res$choice$model_name
    thr  <- unlist(res$summary$threshold)
    lo   <- unlist(res$summary$CI_low)
    hi   <- unlist(res$summary$CI_high)
    tibble(
      Metric                 = metric_names[[metric]],
      `Window size`          = as.integer(gsub("\\D", "", win_nm)),
      `Linear AIC`           = get_aic(res$aic_table, "lm"),
      `Best model AIC`       = if (is.na(best)) NA_real_ else get_aic(res$aic_table, best),
      `Best model`           = if (is.na(best)) "to be decided" else best,
      `Threshold 1 [95% CI]` = fmt_thr(thr, lo, hi, 1),
      `Threshold 2 [95% CI]` = fmt_thr(thr, lo, hi, 2)
    )
  })
}) |>
  mutate(Metric = factor(Metric, levels = metric_names)) |>
  arrange(Metric, `Window size`) |>
  mutate(across(c(`Linear AIC`, `Best model AIC`), \(x) sprintf("%.1f", x)))

latex_table <- kbl(tab, format = "latex", booktabs = TRUE,
                   align = "lcccccc", linesep = "",
                   caption = "Comparison of selected models across window sizes, for each functional diversity metric.") |>
  kable_styling(latex_options = "hold_position", font_size = 9) |>
  collapse_rows(columns = 1, valign = "middle", latex_hline = "major") |>
  landscape()

cat(latex_table)   # à copier dans ton document
# ou pour l'enregistrer :
# writeLines(latex_table, here::here("figures", "article_1", "table_models_windows.tex"))