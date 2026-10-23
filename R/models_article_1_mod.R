# =========================================================================
# Fonctions pour l'analyse de sensibilite (moving window x metrique x modele)
# =========================================================================
#
# Structure generale :
#
#   fit_models()            -> ajuste tous les modeles (SANS bootstrap)
#   compare_models()        -> table AIC / deltaAIC / k (rien d'autre)
#   check_models()          -> diagnostics des residus (hist, fitted vs res, QQ)
#   choose_model()          -> logique "automatique" ou "manuelle" (deltaAIC < 2)
#   refit_with_bootstrap()  -> bootstrap sur UN seul modele (le modele retenu)
#   extract_model_info()    -> threshold + IC, quel que soit le type de modele
#   plot_model()            -> graphique adapte au type de modele
#   analyse_metric()        -> orchestre les fonctions ci-dessus pour 1 metrique
#   analyse_window()        -> boucle analyse_metric() sur les 5 metriques
#   run_sensitivity()       -> lance tout avec barre de progression
#
# Strategie de vitesse : le bootstrap est l'etape la plus lente. Il n'est
# execute que pour le modele retenu (automatique ou impose via
# selected_model). En cas de "manual", aucun bootstrap n'est fait.
#
# Modele stegmented : la recherche du seuil est limitee a 100 points de grille
# (grid.search.max = 100). Sur 4000 lignes, cette grille donne le meme seuil
# et le meme AIC que la recherche complete pour les 15 combinaisons
# fenetre x metrique testees. Le reglage s'applique a l'ajustement rapide ET
# aux replicats du bootstrap. Le nombre de replicats reste celui par defaut
# de chngptm (ci.bootstrap.size = 1000).
#
# Barre de progression : package progressr (une etape par modele ajuste
# + une etape par metrique pour le bootstrap final).
#
# A ADAPTER : les blocs marques "# TODO adapter" dependent de la structure
# exacte de tes objets chngptm (str(model) pour verifier les noms exacts).
# =========================================================================

library(dplyr)
library(tibble)
library(purrr)
library(progressr)
library(chngpt)   # pour chngptm()
# library(mgcv)     # si tu gardes des modeles gam




# -------------------------------------------------------------------------
# 0. Specification des modeles a ajuster
# -------------------------------------------------------------------------
# Chaque element est une fonction(data, response, predictor, bootstrap)
# qui renvoie un modele ajuste. Par defaut (bootstrap = FALSE) : ajustement
# rapide, sans variance. Avec bootstrap = TRUE : bootstrap wildsieve.
#
# make_chngpt_spec(type, ...) : les arguments supplementaires (...) sont
# transmis tels quels a chngptm() (ex : grid.search.max = 100).

make_chngpt_spec <- function(type, ...) {
  extra <- list(...)
  # types pour lesquels chngpt travaille en interne sur -x (hinge.to.upperhinge) :
  # les donnees doivent alors etre triees par ordre decroissant du predicteur,
  # sinon le bootstrap par residus (wild, sieve, wildsieve, awb) apparie
  # les y simules aux mauvais x
  flip_types <- c("hinge", "M01", "M02", "M03", "M04", "M12", "M12c", "M13")
  function(data, response, predictor, bootstrap = FALSE) {
    data <- data[order(data[[predictor]], decreasing = type %in% flip_types), ]
    do.call(chngpt::chngptm, c(list(
      formula.1 = as.formula(paste(response, "~ 1")),
      formula.2 = as.formula(paste("~", predictor)),
      family = "gaussian",
      type = type,
      data = data,
      var.type = if (bootstrap) "bootstrap" else "none",
      bootstrap.type = "awb"
    ), extra))
  }
}

default_model_specs <- function() {
  list(
    lm = function(data, response, predictor, bootstrap = FALSE) {
      lm(as.formula(paste(response, "~", predictor)), data = data)
    },
    step       = make_chngpt_spec("step"),
    segmented  = make_chngpt_spec("segmented"),
    stegmented = make_chngpt_spec("stegmented"),
    M111       = make_chngpt_spec("M111"),
    M02        = make_chngpt_spec("M02"),
    M12        = make_chngpt_spec("M12")
  )
}

# -------------------------------------------------------------------------
# 1. fit_models() -- ajustement rapide, sans bootstrap
# -------------------------------------------------------------------------
# p : fonction de progression (progressr). Par defaut : ne fait rien.

fit_models <- function(data, response, predictor = "median_log_sum",
                       model_specs = default_model_specs(),
                       p = function(...) NULL) {
  
  models <- purrr::imap(model_specs, function(fit_fun, name) {
    p(sprintf("%s | ajustement %s", response, name))
    message(sprintf("[%s] DEBUT %s | %s", format(Sys.time(), "%H:%M:%S"), response, name))
    t0 <- Sys.time()
    out <- tryCatch(
      fit_fun(data, response, predictor, bootstrap = FALSE),
      error = function(e) {
        message(sprintf("Echec ajustement '%s' : %s", name, conditionMessage(e)))
        NULL
      }
    )
    message(sprintf("[%s] FIN   %s | %s : %.1f s", format(Sys.time(), "%H:%M:%S"),
                    response, name, as.numeric(difftime(Sys.time(), t0, units = "secs"))))
    out
  })
  
  models[!vapply(models, is.null, logical(1))]
}

# -------------------------------------------------------------------------
# 1b. refit_with_bootstrap() -- bootstrap sur UN seul modele
# -------------------------------------------------------------------------
refit_with_bootstrap <- function(model_name, data, response, predictor,
                                 model_specs) {
  if (model_name == "lm") {
    return(model_specs[["lm"]](data, response, predictor))
  }
  tryCatch(
    model_specs[[model_name]](data, response, predictor, bootstrap = TRUE),
    error = function(e) {
      message(sprintf("Echec bootstrap '%s' : %s", model_name, conditionMessage(e)))
      NULL
    }
  )
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
# Etapes de progression par metrique : length(model_specs) + 1
# (un ajustement par modele + une etape "bootstrap", meme si sautee).

analyse_metric <- function(data, response, predictor = "median_log_sum",
                           model_specs = default_model_specs(),
                           selected_model = NULL,
                           delta_threshold = 2,
                           run_diagnostics = TRUE,
                           seed = NULL,
                           p = function(...) NULL) {
  
  # 1. Ajustement rapide de tous les modeles (sans bootstrap)
  models      <- fit_models(data, response, predictor, model_specs, p = p)
  aic_table   <- compare_models(models)
  diagnostics <- if (run_diagnostics) check_models(models) else NULL
  choice      <- choose_model(aic_table, models, selected = selected_model,
                              delta_threshold = delta_threshold)
  
  # 2. Bootstrap uniquement sur le modele retenu
  final_model <- NULL
  if (choice$status != "manual") {
    p(sprintf("%s | bootstrap %s", response, choice$model_name))
    message(sprintf("[%s] DEBUT bootstrap %s | %s", format(Sys.time(), "%H:%M:%S"),
                    response, choice$model_name))
    if (!is.null(seed)) set.seed(seed)
    final_model <- refit_with_bootstrap(choice$model_name, data, response,
                                        predictor, model_specs)
    message(sprintf("[%s] FIN   bootstrap %s | %s", format(Sys.time(), "%H:%M:%S"),
                    response, choice$model_name))
    # si le bootstrap echoue, on garde le modele sans bootstrap
    if (is.null(final_model)) final_model <- choice$model
  } else {
    p(sprintf("%s | manuel : pas de bootstrap", response))
  }
  
  info <- if (!is.null(final_model)) {
    extract_model_info(final_model, choice$model_name)
  } else {
    tibble::tibble(model = NA_character_, threshold = NA_real_,
                   CI_low = NA_real_, CI_high = NA_real_)
  }
  
  plot <- if (!is.null(final_model)) {
    plot_model(final_model, data, response, predictor, choice$model_name)
  } else {
    NULL
  }
  
  list(
    models      = models,        # ajustes sans bootstrap
    aic_table   = aic_table,
    diagnostics = diagnostics,
    choice      = choice,
    final_model = final_model,   # modele retenu, avec bootstrap
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
                           run_diagnostics = TRUE,
                           seed = NULL,
                           p = function(...) NULL) {
  
  purrr::map(rlang::set_names(responses), function(resp) {
    analyse_metric(
      data            = data,
      response        = resp,
      predictor       = predictor,
      model_specs     = model_specs,
      selected_model  = selected_models[[resp]],
      delta_threshold = delta_threshold,
      run_diagnostics = run_diagnostics,
      seed            = seed,
      p               = p
    )
  })
}

# -------------------------------------------------------------------------
# 9. run_sensitivity() -- lance plusieurs fenetres avec barre de progression
# -------------------------------------------------------------------------
# windows : liste nommee de data.frames, ex:
#   list(window30 = sum_df_30_25draw, window35 = sum_df_35_25draw, ...)

run_sensitivity <- function(windows,
                            responses = c("Richness", "regularity", "functional_dispersion",
                                          "originality", "dimensionality"),
                            predictor = "median_log_sum",
                            model_specs = default_model_specs(),
                            selected_models = list(),
                            delta_threshold = 2,
                            run_diagnostics = TRUE,
                            seed = NULL) {
  
  n_steps <- length(windows) * length(responses) * (length(model_specs) + 1)
  
  progressr::with_progress({
    p <- progressr::progressor(steps = n_steps)
    purrr::imap(windows, function(df, nm) {
      analyse_window(
        df,
        responses       = responses,
        predictor       = predictor,
        model_specs     = model_specs,
        selected_models = selected_models,
        delta_threshold = delta_threshold,
        run_diagnostics = run_diagnostics,
        seed            = seed,
        p               = p
      )
    })
  })
}

# =========================================================================
# Exemple d'utilisation (a lancer dans la console, pas dans ce fichier)
# =========================================================================
#
 progressr::handlers(global = TRUE)
 progressr::handlers("cli")   # affiche aussi le message ; sinon "txtprogressbar"
#
 results <- run_sensitivity(
   windows = list(window30 = sum_df_30_25draw_q95 ,
                 window35 = sum_df_35_25draw_q95,
                  window40 = sum_df_40_25draw_q95),
   seed = 1
 )
#
# results$window30$Richness$choice$status   # "automatic", "manual" ou "manual_forced"
# results$window30$Richness$summary
#
# # cas manuel : tu regardes les candidats et tu relances avec ton choix
# results$window30$Richness$choice$candidates
# results$window30$Richness <- analyse_metric(
#   sum_df_30_25draw, "Richness", selected_model = "segmented", seed = 1
# )
 
 