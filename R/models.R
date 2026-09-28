fit_feature <- function(y, md, cfg, sensitivity = FALSE) {
  d <- md
  d$response <- as.numeric(y)
  warnings <- character()
  if (stats::sd(d$response) < 1e-10) return(list(error = "No response variation", warnings = warnings))
  fit <- tryCatch(withCallingHandlers(
    nlme::lme(fixed = model_formula(cfg, sensitivity), random = ~ 1 | participant_id,
      data = d, method = "REML", na.action = stats::na.fail,
      correlation = if (cfg$correlation == "car1") nlme::corCAR1(form = ~ time_scaled | participant_id) else NULL,
      control = nlme::lmeControl(maxIter = 100, msMaxIter = 200, returnObject = FALSE)),
    warning = function(w) { warnings <<- c(warnings, conditionMessage(w)); invokeRestart("muffleWarning") }),
    error = function(e) e)
  if (inherits(fit, "error")) return(list(error = conditionMessage(fit), warnings = warnings))
  tt <- summary(fit)$tTable
  if (any(!is.finite(tt[,c("Value","Std.Error","DF","p-value")])))
    return(list(error = "Nonfinite model statistics", warnings = warnings))
  random_sd <- as.numeric(nlme::VarCorr(fit)[1,"StdDev"])
  if (random_sd / fit$sigma < 1e-4) warnings <- c(warnings, "Random intercept variance near boundary")
  list(table = tt, residuals = as.numeric(stats::residuals(fit, type = "normalized")),
    fitted = as.numeric(stats::fitted(fit)), sigma = fit$sigma, random_sd = random_sd,
    rho = if (cfg$correlation == "car1") as.numeric(stats::coef(fit$modelStruct$corStruct, unconstrained = FALSE)) else NA_real_,
    warnings = unique(warnings))
}

fit_modality <- function(data, md, cfg, label, sensitivity = FALSE) {
  coefficients <- list(); diagnostics <- list(); residuals <- list(); fitted <- list()
  spec <- if (sensitivity) "sensitivity" else "primary"
  for (feature in data$keep) {
    result <- fit_feature(data$clr[, feature], md, cfg, sensitivity)
    status <- if (is.null(result$error)) "ok" else "failed"
    diagnostics[[feature]] <- data.frame(modality = label, feature = feature, specification = spec,
      status = status, message = paste(c(result$error,result$warnings), collapse = "; "),
      n_samples = nrow(md), n_participants = nlevels(md$participant_id),
      residual_sd = if (status == "ok") result$sigma else NA_real_,
      random_intercept_sd = if (status == "ok") result$random_sd else NA_real_,
      car1_rho_per_time_unit = if (status == "ok") result$rho else NA_real_)
    if (status != "ok") next
    tt <- result$table
    terms <- rownames(tt)
    primary <- if (is.null(cfg$group_reference)) grepl("^phase",terms) else grepl(":",terms,fixed=TRUE)
    se <- tt[,"Std.Error"]; est <- tt[,"Value"]; crit <- stats::qt(0.975, tt[,"DF"])
    coefficients[[feature]] <- data.frame(modality = label, feature = feature, specification = spec,
      term = terms, primary_contrast = primary, estimate = est, standard_error = se,
      df = tt[,"DF"], ci_lower = est - crit * se, ci_upper = est + crit * se,
      p_value = tt[,"p-value"], q_value = NA_real_, stringsAsFactors = FALSE, row.names = NULL)
    residuals[[feature]] <- result$residuals
    fitted[[feature]] <- result$fitted
  }
  list(coefficients = if (length(coefficients)) do.call(rbind, coefficients) else NULL,
    diagnostics = do.call(rbind, diagnostics), residuals = residuals, fitted = fitted)
}

fit_study <- function(study) {
  cfg <- study$config
  models <- list()
  for (mod in c("bacteria", if (!is.null(study$phages)) "phages")) {
    models[[mod]] <- fit_modality(study[[mod]], study$metadata, cfg, mod)
    if (length(cfg$sensitivity_covariates)) models[[paste0(mod,"_sensitivity")]] <-
      fit_modality(study[[mod]], study$metadata, cfg, mod, TRUE)
  }
  coefficients <- do.call(rbind, lapply(models, `[[`, "coefficients"))
  diagnostics <- do.call(rbind, lapply(models, `[[`, "diagnostics"))
  if (is.null(coefficients)) cb_stop("Every mixed model failed; no report was generated. Check study design and response variation.")
  # Correct all planned feature-by-phase contrasts jointly across modalities,
  # separately by specification. Failed fits remain in the family size.
  for (spec in unique(coefficients$specification)) {
    idx <- coefficients$specification == spec & coefficients$primary_contrast
    planned <- sum(diagnostics$specification == spec) * (length(cfg$phase_levels)-1L)
    coefficients$q_value[idx] <- stats::p.adjust(coefficients$p_value[idx], method="BH", n = planned)
  }
  rownames(coefficients) <- NULL; rownames(diagnostics) <- NULL
  list(models = models, coefficients = coefficients, diagnostics = diagnostics)
}

with_cb_seed <- function(seed, fun) {
  had_seed <- exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
  if (had_seed) old_seed <- get(".Random.seed", envir = .GlobalEnv)
  on.exit(if (had_seed) assign(".Random.seed", old_seed, envir = .GlobalEnv) else
    if (exists(".Random.seed", envir = .GlobalEnv, inherits=FALSE)) rm(".Random.seed", envir = .GlobalEnv))
  set.seed(seed)
  fun()
}

host_associations <- function(study, fits) {
  empty <- data.frame(phage_id=character(), host_id=character(), confidence=character(), evidence=character(),
    status=character(), n_participants=integer(), mean_within_person_rho=numeric(),
    ci_lower=numeric(), ci_upper=numeric(), included_primary=logical(), included_high_only=logical())
  if (is.null(study$links)) return(empty)
  with_cb_seed(study$config$seed, function() {
    rows <- lapply(seq_len(nrow(study$links)), function(i) {
      link <- study$links[i,]
      p <- fits$models$phages$residuals[[link$phage_id]]
      b <- fits$models$bacteria$residuals[[link$host_id]]
      rhos <- numeric()
      if (!is.null(p) && !is.null(b)) for (idx in split(seq_len(nrow(study$metadata)), study$metadata$participant_id)) {
        if (length(idx) >= 4 && stats::sd(p[idx]) > 1e-10 && stats::sd(b[idx]) > 1e-10)
          rhos <- c(rhos, stats::cor(p[idx], b[idx], method="spearman"))
      }
      rhos <- rhos[is.finite(rhos)]
      ok <- length(rhos) >= 6
      ci <- if (ok) stats::quantile(replicate(study$config$bootstrap_replicates,
        mean(sample(rhos, length(rhos), replace = TRUE))), c(0.025,0.975), names=FALSE) else c(NA_real_,NA_real_)
      data.frame(phage_id=link$phage_id, host_id=link$host_id, confidence=link$confidence,
        evidence=link$evidence, status=if (is.null(p) || is.null(b)) "model_unavailable_or_filtered" else
          if (ok) "exploratory" else "insufficient_participants",
        n_participants=length(rhos), mean_within_person_rho=if (ok) mean(rhos) else NA_real_,
        ci_lower=ci[1], ci_upper=ci[2], included_primary=link$confidence %in% study$config$host_confidence,
        included_high_only=link$confidence == "high", stringsAsFactors=FALSE)
    })
    do.call(rbind, rows)
  })
}
