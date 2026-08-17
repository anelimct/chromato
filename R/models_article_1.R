# =========================================================================
# Fonctions pour l'analyse de sensibilite (moving window x metrique x modele)
# =========================================================================
#
# Structure generale :
#
#   fit_models()        -> ajuste tous les modeles pour une metrique
#   compare_models()     -> table AIC / deltaAIC / k (rien d'autre)
#   check_models()        -> diagnostics des residus (hist, fitted vs res, QQ)
#   choose_model()        -> logique "automatique" ou "manuelle" (deltaAIC < 2)
#   extract_model_info()  -> threshold + IC, quel que soit le type de modele
#   plot_model()          -> graphique adapte au type de modele
#   analyse_metric()      -> orchestre les 5 fonctions ci-dessus pour 1 metrique
#   analyse_window()      -> boucle analyse_metric() sur les 5 metriques
#
# A ADAPTER : les blocs marques "# TODO adapter" dependent de la structure
# exacte de tes objets chngptm (str(model) pour verifier les noms exacts).
# =========================================================================

library(dplyr)
library(tibble)
library(purrr)
# library(chngpt)   # pour chngptm()
# library(mgcv)     # si tu gardes des modeles gam



# -------------------------------------------------------------------------
# 0. Specification des modeles a ajuster
# -------------------------------------------------------------------------
# Chaque element est une fonction(data, response, predictor) qui renvoie
# un modele ajuste (ou NULL / un objet try-error si l'ajustement echoue).
# Adapte les appels chngptm() ci-dessous a ta formule/famille exacte.

default_model_specs <- function() {
  list(
    lm = function(data, response, predictor) {
      f <- as.formula(paste(response, "~", predictor))
      lm(f, data = data)
    },
    
    step = function(data, response, predictor) {
      # TODO adapter : formule / family
      chngpt::chngptm(
        formula.1 = as.formula(paste(response, "~ 1")),
        formula.2 = as.formula(paste("~", predictor)),
        family = "gaussian",
        type = "step",
        data = data, var.type = "bootstrap", bootsrap.type = "wildsieve" 
      )
    },
    
    segmented = function(data, response, predictor) {
      chngpt::chngptm(
        formula.1 = as.formula(paste(response, "~ 1")),
        formula.2 = as.formula(paste("~", predictor)),
        family = "gaussian",
        type = "segmented",
        data = data, var.type = "bootstrap", bootsrap.type = "wildsieve"
      )
    },
    
    stegmented = function(data, response, predictor) {
      chngpt::chngptm(
        formula.1 = as.formula(paste(response, "~ 1")),
        formula.2 = as.formula(paste("~", predictor)),
        family = "gaussian",
        type = "stegmented",
        data = data, var.type = "bootstrap", bootsrap.type = "wildsieve"
      )
    },
    
    M111 = function(data, response, predictor) {
      chngpt::chngptm(
        formula.1 = as.formula(paste(response, "~ 1")),
        formula.2 = as.formula(paste("~", predictor)),
        family = "gaussian",
        type = "M111",
        data = data, var.type = "bootstrap", bootsrap.type = "wildsieve"
      )
    },
    
    M02 = function(data, response, predictor) {
      chngpt::chngptm(
        formula.1 = as.formula(paste(response, "~ 1")),
        formula.2 = as.formula(paste("~", predictor)),
        family = "gaussian",
        type = "M02",
        data = data, var.type = "bootstrap", bootsrap.type = "wildsieve"
      )
    },
    
    M12 = function(data, response, predictor) {
      chngpt::chngptm(
        formula.1 = as.formula(paste(response, "~ 1")),
        formula.2 = as.formula(paste("~", predictor)),
        family = "gaussian",
        type = "M12",
        data = data, var.type = "bootstrap", bootsrap.type = "wildsieve"
      )
    }

  )
}

# -------------------------------------------------------------------------
# 1. fit_models()
# -------------------------------------------------------------------------
fit_models <- function(data, response, predictor = "median_log_sum",
                       model_specs = default_model_specs()) {
  
  models <- purrr::imap(model_specs, function(fit_fun, name) {
    out <- tryCatch(
      fit_fun(data, response, predictor),
      error = function(e) {
        message(sprintf("Echec ajustement '%s' : %s", name, conditionMessage(e)))
        NULL
      }
    )
    out
  })
  
  models[!vapply(models, is.null, logical(1))]
}

# -------------------------------------------------------------------------
# 2. compare_models()  -- AIC / deltaAIC / k UNIQUEMENT
# -------------------------------------------------------------------------
safe_AIC <- function(model) {
  tryCatch(AIC(model), error = function(e) NA_real_)
}

safe_k <- function(model) {
  tryCatch({
    cf <- coef(model)
    if (is.list(cf)) length(unlist(cf)) else length(cf)
  }, error = function(e) NA_integer_)
}

compare_models <- function(models) {
  tibble::tibble(
    model = names(models),
    AIC   = vapply(models, safe_AIC, numeric(1)),
    k     = vapply(models, safe_k, integer(1))
  ) %>%
    dplyr::arrange(AIC) %>%
    dplyr::mutate(deltaAIC = AIC - min(AIC, na.rm = TRUE)) %>%
    dplyr::relocate(model, AIC, deltaAIC, k)
}

# -------------------------------------------------------------------------
# 3. check_models() -- diagnostics des residus, separe de compare_models()
# -------------------------------------------------------------------------
get_residuals <- function(model) {
  if (inherits(model, "lm")) {
    return(list(fitted = model$fitted.values, resid = rstandard(model),
                standardized = TRUE))
  }
  if (inherits(model, "chngptm")) {
    # structure reelle : les residus/fitted du modele final sont dans best.fit
    fitted_vals <- tryCatch(model$best.fit$fitted.values, error = function(e) NULL)
    resid_vals  <- tryCatch(model$best.fit$residuals,     error = function(e) NULL)
    
    # filet de securite si best.fit n'existe pas pour un type de modele donne
    if (is.null(resid_vals)) {
      resid_vals <- tryCatch(residuals(model), error = function(e) NULL)
    }
    if (is.null(fitted_vals)) {
      fitted_vals <- tryCatch(fitted(model), error = function(e) NULL)
    }
    
    return(list(fitted = fitted_vals, resid = resid_vals, standardized = FALSE))
  }
  list(fitted = NULL, resid = NULL, standardized = FALSE)
}

check_model <- function(model, model_name = "") {
  rr <- get_residuals(model)
  if (is.null(rr$resid) || is.null(rr$fitted)) {
    message(sprintf("Pas de residus extractibles pour '%s'", model_name))
    return(invisible(NULL))
  }
  
  op <- par(mfrow = c(1, 3))
  on.exit(par(op))
  
  ylab_txt <- if (isTRUE(rr$standardized)) "Standardized residuals" else "Residuals"
  
  hist(rr$resid, main = paste(model_name, "- residuals"), col = "#889cca",
       xlab = ylab_txt)
  
  plot(rr$fitted, rr$resid,
       main = paste(model_name, "-", ylab_txt, "vs Fitted"),
       xlab = "Fitted values", ylab = ylab_txt,
       col = "#889cca")
  abline(h = 0, col = "grey")
  
  qqnorm(rr$resid, main = paste(model_name, "- QQ-plot"), col = "#889cca")
  qqline(rr$resid, col = "grey")
  
  invisible(list(
    shapiro_p = tryCatch(shapiro.test(rr$resid)$p.value, error = function(e) NA_real_)
  ))
}

check_models <- function(models) {
  purrr::imap(models, function(m, nm) check_model(m, nm))
}


# -------------------------------------------------------------------------
# 4. choose_model() -- logique automatique / manuelle
# -------------------------------------------------------------------------
choose_model <- function(aic_table, models, selected = NULL, delta_threshold = 2) {
  
  if (!is.null(selected)) {
    return(list(
      status = "manual_forced",
      model_name = selected,
      model = models[[selected]],
      candidates = aic_table
    ))
  }
  
  best_name <- aic_table$model[1]
  ambiguous <- aic_table %>% dplyr::filter(deltaAIC < delta_threshold)
  
  if (nrow(ambiguous) > 1) {
    return(list(
      status = "manual",           # toi seul tranches
      model_name = NA_character_,
      model = NULL,
      candidates = ambiguous
    ))
  }
  
  list(
    status = "automatic",
    model_name = best_name,
    model = models[[best_name]],
    candidates = aic_table[1, ]
  )
}

# -------------------------------------------------------------------------
# 5. extract_model_info() -- threshold + IC, masque les differences de classe
# -------------------------------------------------------------------------
extract_model_info <- function(model, model_name = "") {
  
  if (is.null(model) || inherits(model, "lm")) {
    return(tibble::tibble(
      model = model_name, threshold = NA_real_,
      CI_low = NA_real_, CI_high = NA_real_
    ))
  }
  
  if (inherits(model, "chngptm")) {
    
    thr   <- tryCatch(as.numeric(model$chngpt), error = function(e) NA_real_)
    n_thr <- length(thr)
    
    ci_low  <- rep(NA_real_, n_thr)
    ci_high <- rep(NA_real_, n_thr)
    found   <- FALSE
    
    # ----- Methode 1 : model$vcov$perc (bootstrap percentiles) -----
    perc_mat <- tryCatch(model$vcov$perc, error = function(e) NULL)
    
    if (!is.null(perc_mat) && is.matrix(perc_mat)) {
      col_names <- colnames(perc_mat)
      chngpt_cols <- grep("^chngpt(\\.[0-9]+)?$", col_names, value = TRUE)
      chngpt_cols <- chngpt_cols[order(as.numeric(gsub("^chngpt\\.?", "0", chngpt_cols)))]
      
      if (length(chngpt_cols) == n_thr &&
          all(c("2.5%", "97.5%") %in% rownames(perc_mat))) {
        ci_low  <- unname(perc_mat["2.5%",  chngpt_cols])
        ci_high <- unname(perc_mat["97.5%", chngpt_cols])
        found <- TRUE
      }
    }
    
    # ----- Methode 2 (fallback) : summary(model) -----
    if (!found) {
      s <- tryCatch(summary(model), error = function(e) NULL)
      if (!is.null(s)) {
        candidate_names <- c("chngpt", "coef.chngpt", "chngpt.coef", "threshold")
        for (nm in candidate_names) {
          val <- tryCatch(s[[nm]], error = function(e) NULL)
          if (!is.null(val)) {
            val_mat <- if (is.matrix(val)) val else matrix(val, nrow = 1)
            if (nrow(val_mat) == n_thr && ncol(val_mat) >= 4) {
              cn <- colnames(val_mat)
              lo_col <- if (!is.null(cn) && "lower" %in% cn) which(cn == "lower") else 3
              hi_col <- if (!is.null(cn) && "upper" %in% cn) which(cn == "upper") else 4
              ci_low  <- val_mat[, lo_col]
              ci_high <- val_mat[, hi_col]
              found <- TRUE
              break
            }
          }
        }
      }
    }
    
    if (!found) {
      message(sprintf(
        "extract_model_info('%s') : IC du threshold introuvable (ni vcov$perc, ni summary())",
        model_name
      ))
    }
    
    return(tibble::tibble(
      model     = model_name,
      threshold = list(thr),
      CI_low    = list(ci_low),
      CI_high   = list(ci_high)
    ))
  }
  
  tibble::tibble(
    model = model_name, threshold = NA_real_,
    CI_low = NA_real_, CI_high = NA_real_
  )
}

# -------------------------------------------------------------------------
# 6. plot_model() -- graphique adapte au type de modele
# -------------------------------------------------------------------------
plot_model <- function(model, data, response, predictor = "median_log_sum", model_name = "") {
  
  if (inherits(model, "lm")) {
    plot(data[[predictor]], data[[response]],
         main = model_name, xlab = predictor, ylab = response, col = "#889cca")
    abline(model, col = "darkred", lwd = 2)
    return(invisible(NULL))
  }
  
  if (inherits(model, "chngptm")) {
    # le package chngpt fournit generalement une methode plot.chngptm()
    tryCatch(
      plot(model, main = model_name),
      error = function(e) message("Pas de methode plot() disponible pour ce modele")
    )
    return(invisible(NULL))
  }
  
  message(sprintf("plot_model() : classe non geree pour '%s'", model_name))
}

# -------------------------------------------------------------------------
# 7. analyse_metric() -- orchestre tout pour UNE metrique
# -------------------------------------------------------------------------
analyse_metric <- function(data, response, predictor = "median_log_sum",
                           model_specs = default_model_specs(),
                           selected_model = NULL,
                           delta_threshold = 2,
                           run_diagnostics = TRUE) {
  
  models     <- fit_models(data, response, predictor, model_specs)
  aic_table  <- compare_models(models)
  diagnostics <- if (run_diagnostics) check_models(models) else NULL
  choice     <- choose_model(aic_table, models, selected = selected_model,
                             delta_threshold = delta_threshold)
  
  info <- if (choice$status != "manual") {
    extract_model_info(choice$model, choice$model_name)
  } else {
    tibble::tibble(model = NA_character_, threshold = NA_real_,
                   CI_low = NA_real_, CI_high = NA_real_)
  }
  
  plot <- if (choice$status != "manual") {
    plot_model(choice$model, data, response, predictor, choice$model_name)
  } else {
    NULL
  }
  
  list(
    models      = models,
    aic_table   = aic_table,
    diagnostics = diagnostics,
    choice      = choice,
    summary     = info,
    plot        = plot
  )
}

# -------------------------------------------------------------------------
# 8. analyse_window() -- boucle analyse_metric() sur toutes les metriques
# -------------------------------------------------------------------------
analyse_window <- function(data, responses = c("Richness", "regularity", "functional_dispersion",
                                               "originality", "dimensionality"),
                           predictor = "median_log_sum",
                           model_specs = default_model_specs(),
                           selected_models = list(),   # ex: list(Richness = "segmented")
                           delta_threshold = 2,
                           run_diagnostics = TRUE) {
  
  purrr::map(rlang::set_names(responses), function(resp) {
    analyse_metric(
      data            = data,
      response        = resp,
      predictor       = predictor,
      model_specs     = model_specs,
      selected_model  = selected_models[[resp]],
      delta_threshold = delta_threshold,
      run_diagnostics = run_diagnostics
    )
  })
}

# =========================================================================
# Exemple d'utilisation
# =========================================================================
#
 results <- list(
   window15 = analyse_window(sum_df_15),
   window20 = analyse_window(sum_df_20),
   window25 = analyse_window(sum_df_25),
   window30 = analyse_window(sum_df_30),
   window35 = analyse_window(sum_df_35),
   window40 = analyse_window(sum_df_40)
 )
#
# 
# results$window25$Richness$choice$status   # "automatic" ou "manual"
# results$window25$Richness$summary
#
# # cas manuel : tu regardes les candidats et tu relances avec ton choix
# results$window25$Richness$choice$candidates
# results$window25$Richness <- analyse_metric(
#   sum_df_25, "Richness", selected_model = "segmented"
# )