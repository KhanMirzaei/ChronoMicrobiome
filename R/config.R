cb_stop <- function(...) stop(..., call. = FALSE)

cb_file <- function(...) {
  local_root <- getOption("chronomicrobiome.root")
  if (!is.null(local_root)) return(file.path(local_root, "inst", ...))
  system.file(..., package = "chronomicrobiome", mustWork = TRUE)
}

read_config <- function(path) {
  path <- normalizePath(path, mustWork = TRUE)
  cfg <- yaml::read_yaml(path, eval.expr = FALSE)
  if (!is.list(cfg)) cb_stop("Configuration must be a YAML mapping.")
  defaults <- list(title = "Chrono\u00b5Biome study", synthetic = FALSE,
    abundance_scale = "counts", phase_levels = c("baseline", "treatment", "followup"),
    group_reference = NULL, covariates = character(), sensitivity_covariates = character(),
    min_prevalence = 0.1, pseudocount = 1e-6, correlation = "independent",
    time_unit_days = 7, recovery_threshold = 0.15, recovery_consecutive = 2L,
    min_baseline_samples = 2L, detection_threshold = 0,
    host_confidence = c("high", "medium"), bootstrap_replicates = 500L,
    seed = 20260928L, report_engine = "auto")
  allowed <- c(names(defaults), "metadata", "bacteria", "phages", "host_links", "output")
  unknown <- setdiff(names(cfg), allowed)
  if (length(unknown)) cb_stop("Unknown configuration keys: ", paste(unknown, collapse = ", "))
  for (nm in names(defaults)) if (is.null(cfg[[nm]])) cfg[[nm]] <- defaults[[nm]]
  for (nm in c("phase_levels", "covariates", "sensitivity_covariates", "host_confidence"))
    cfg[[nm]] <- as.character(unlist(cfg[[nm]], use.names = FALSE))
  if (!is.logical(cfg$synthetic) || length(cfg$synthetic) != 1 || is.na(cfg$synthetic))
    cb_stop("synthetic must be true or false.")
  for (nm in c("title", "abundance_scale", "correlation", "report_engine"))
    if (!is.character(cfg[[nm]]) || length(cfg[[nm]]) != 1 || is.na(cfg[[nm]])) cb_stop(nm, " must be a string.")
  if (!cfg$abundance_scale %in% c("counts", "relative")) cb_stop("abundance_scale must be counts or relative.")
  if (!cfg$correlation %in% c("independent", "car1")) cb_stop("correlation must be independent or car1.")
  if (!cfg$report_engine %in% c("auto", "quarto", "builtin")) cb_stop("report_engine must be auto, quarto, or builtin.")
  if (length(cfg$phase_levels) < 2 || anyNA(cfg$phase_levels) || any(!nzchar(cfg$phase_levels)) ||
      anyDuplicated(cfg$phase_levels) || cfg$phase_levels[1] != "baseline")
    cb_stop("phase_levels must start with baseline and contain at least two unique levels.")
  bounds <- list(min_prevalence = c(0,1), pseudocount = c(1e-12,0.01),
    time_unit_days = c(1e-6,1e6), recovery_threshold = c(0,1), recovery_consecutive = c(1,100),
    min_baseline_samples = c(1,100), detection_threshold = c(0,1),
    bootstrap_replicates = c(100,100000), seed = c(0,.Machine$integer.max))
  for (nm in names(bounds)) {
    v <- cfg[[nm]]
    if (!is.numeric(v) || length(v) != 1 || !is.finite(v) || v < bounds[[nm]][1] || v > bounds[[nm]][2])
      cb_stop("Invalid numeric configuration: ", nm)
  }
  for (nm in c("recovery_consecutive", "min_baseline_samples", "bootstrap_replicates", "seed"))
    if (cfg[[nm]] != as.integer(cfg[[nm]])) cb_stop(nm, " must be an integer.")
  if (!length(cfg$host_confidence) || any(!cfg$host_confidence %in% c("high","medium","low")))
    cb_stop("host_confidence must contain high, medium, or low.")
  reserved <- c("sample_id", "participant_id", "time", "phase", "group", "response", "time_scaled")
  cv <- c(cfg$covariates, cfg$sensitivity_covariates)
  if (anyDuplicated(cv) || any(!grepl("^[A-Za-z][A-Za-z0-9_]*$", cv)) || any(cv %in% reserved))
    cb_stop("Covariates must be unique simple column names, excluding core metadata columns.")
  if (!is.null(cfg$group_reference) && (!is.character(cfg$group_reference) || length(cfg$group_reference) != 1 || !nzchar(cfg$group_reference)))
    cb_stop("group_reference must be a nonempty string or null.")
  for (nm in c("metadata", "bacteria")) if (is.null(cfg[[nm]])) cb_stop("Missing configuration: ", nm)
  base <- dirname(path)
  resolve <- function(x) {
    if (is.null(x)) return(NULL)
    if (!is.character(x) || length(x) != 1 || !nzchar(x)) cb_stop("File paths must be nonempty strings.")
    x <- path.expand(x)
    if (!grepl("^(/|[A-Za-z]:[/\\\\])", x)) x <- file.path(base, x)
    x
  }
  for (nm in c("metadata","bacteria","phages","host_links","output")) cfg[[nm]] <- resolve(cfg[[nm]])
  if (is.null(cfg$output)) cfg$output <- file.path(base, "results")
  if (!is.null(cfg$host_links) && is.null(cfg$phages)) cb_stop("host_links requires a phages table.")
  cfg$config_file <- path
  cfg
}

read_table <- function(path) {
  if (!file.exists(path)) cb_stop("Input file does not exist: ", path)
  sep <- if (grepl("\\.csv$", path, ignore.case = TRUE)) "," else "\t"
  x <- utils::read.table(path, header = TRUE, sep = sep, quote = "\"", comment.char = "",
    check.names = FALSE, stringsAsFactors = FALSE, na.strings = c("NA", ""))
  if (!nrow(x) || anyDuplicated(names(x)) || any(!nzchar(names(x)))) cb_stop("Empty table or duplicate/empty columns: ", path)
  x
}

model_formula <- function(cfg, sensitivity = FALSE) {
  rhs <- if (is.null(cfg$group_reference)) "phase" else "group * phase"
  covars <- c(cfg$covariates, if (sensitivity) cfg$sensitivity_covariates)
  stats::as.formula(paste("response ~", paste(c(rhs, covars), collapse = " + ")))
}

read_abundance <- function(path, metadata, cfg, label) {
  x <- read_table(path)
  if (names(x)[1] != "sample_id" || ncol(x) < 3) cb_stop(label, ": first column must be sample_id, followed by at least two features.")
  if (anyNA(x$sample_id) || anyDuplicated(x$sample_id)) cb_stop(label, ": missing or duplicate sample IDs.")
  if (!setequal(x$sample_id, metadata$sample_id)) cb_stop(label, ": sample IDs must exactly match metadata. No samples are silently dropped.")
  x <- x[match(metadata$sample_id, x$sample_id), -1, drop = FALSE]
  if (!all(vapply(x, is.numeric, logical(1)))) cb_stop(label, ": abundances must be numeric.")
  m <- as.matrix(x)
  if (any(!is.finite(m)) || any(m < 0)) cb_stop(label, ": abundances must be finite, nonnegative, and nonmissing.")
  if (cfg$abundance_scale == "counts" && any(abs(m - round(m)) > 1e-8)) cb_stop(label, ": count inputs must be integers; use relative for proportions.")
  totals <- rowSums(m)
  if (any(totals <= 0)) cb_stop(label, ": zero-total samples cannot be normalized.")
  if (cfg$abundance_scale == "relative" && any(abs(totals - 1) > 0.01))
    cb_stop(label, ": relative inputs must sum to 1 per sample (within 0.01), not percentages.")
  rel <- m / totals
  nonzero <- colSums(m) > 0
  if (sum(nonzero) < 2) cb_stop(label, ": at least two nonzero features are required for CLR analysis.")
  rel <- rel[, nonzero, drop = FALSE]
  rownames(rel) <- metadata$sample_id
  logrel <- log(rel + cfg$pseudocount)
  clr <- logrel - rowMeans(logrel)
  prevalence <- colMeans(rel > cfg$detection_threshold)
  keep <- prevalence >= cfg$min_prevalence
  if (!any(keep)) cb_stop(label, ": no features pass min_prevalence.")
  audit <- data.frame(modality = label, feature = colnames(m), prevalence = colMeans(m / totals > cfg$detection_threshold),
    status = ifelse(!nonzero, "all_zero", ifelse(colMeans(m / totals > cfg$detection_threshold) >= cfg$min_prevalence, "modeled", "below_prevalence")))
  list(relative = rel, clr = clr, keep = colnames(rel)[keep], audit = audit, totals = totals)
}

validate_study <- function(config) {
  cfg <- read_config(config)
  md <- read_table(cfg$metadata)
  required <- c("sample_id", "participant_id", "time", "phase", cfg$covariates, cfg$sensitivity_covariates,
    if (!is.null(cfg$group_reference)) "group")
  missing <- setdiff(required, names(md))
  if (length(missing)) cb_stop("Missing metadata columns: ", paste(missing, collapse = ", "))
  if (anyNA(md[, required, drop = FALSE])) cb_stop("Missing values in required metadata. Resolve them explicitly before analysis.")
  if (anyDuplicated(md$sample_id)) cb_stop("Duplicate metadata sample_id.")
  for (nm in c("sample_id", "participant_id", "phase")) {
    md[[nm]] <- as.character(md[[nm]])
    if (any(!nzchar(trimws(md[[nm]])))) cb_stop("Empty metadata identifiers or phase.")
  }
  if (!is.numeric(md$time) || any(!is.finite(md$time))) cb_stop("time must be finite numeric days relative to the study index event.")
  if (anyDuplicated(md[, c("participant_id","time")])) cb_stop("Duplicate participant/time observations. Resolve technical replicates explicitly.")
  if (!setequal(unique(md$phase), cfg$phase_levels)) cb_stop("Observed phase labels must exactly match phase_levels.")
  md$phase <- factor(md$phase, levels = cfg$phase_levels)
  md$participant_id <- factor(md$participant_id)
  md <- md[order(md$participant_id, md$time), , drop = FALSE]
  rownames(md) <- NULL
  if (nlevels(md$participant_id) < 6) cb_stop("At least 6 participants are required by this exploratory workflow; this is not a power guarantee.")
  if (any(table(md$participant_id) < 2)) cb_stop("Each participant must have at least two samples.")
  if (any(vapply(split(as.integer(md$phase), md$participant_id), function(x) any(diff(x) < 0), logical(1))))
    cb_stop("Phases must follow phase_levels within each participant. Recurrent exposure episodes are not supported in v0.1.")
  if (!is.null(cfg$group_reference)) {
    md$group <- as.character(md$group)
    groups <- unique(md$group)
    if (length(groups) != 2 || !cfg$group_reference %in% groups) cb_stop("Two study groups are required, including group_reference.")
    if (any(vapply(split(md$group, md$participant_id), function(x) length(unique(x)) != 1, logical(1))))
      cb_stop("group must remain constant within participant. Use phase for time-varying exposure windows.")
    md$group <- factor(md$group, levels = c(cfg$group_reference, setdiff(groups,cfg$group_reference)))
    participants <- unique(md[, c("participant_id","group")])
    if (any(table(participants$group) < 3)) cb_stop("At least 3 participants per group are required; this is not a power guarantee.")
    if (any(table(md$group, md$phase) == 0)) cb_stop("Each group must be observed in every phase.")
  } else {
    md$group <- factor("single_cohort")
  }
  for (nm in c(cfg$covariates,cfg$sensitivity_covariates)) {
    if (is.character(md[[nm]])) md[[nm]] <- factor(md[[nm]])
    if (is.numeric(md[[nm]]) && any(!is.finite(md[[nm]]))) cb_stop("Nonfinite covariate: ", nm)
    if (length(unique(md[[nm]])) < 2) cb_stop("Covariate has no variation: ", nm)
  }
  md$time_scaled <- md$time / cfg$time_unit_days
  for (sens in c(FALSE, if (length(cfg$sensitivity_covariates)) TRUE)) {
    mm <- stats::model.matrix(stats::delete.response(stats::terms(model_formula(cfg, sens))), md)
    if (qr(mm)$rank < ncol(mm)) cb_stop("Model is rank deficient. Check collinear covariates and phase/group coverage.")
    if (nrow(mm) <= ncol(mm) + nlevels(md$participant_id)) cb_stop("Insufficient observations for the requested model.")
  }
  notes <- c("Exploratory CLR mixed models estimate relative associations, not absolute abundance or causal drug effects.",
    "Detection depends on sequencing depth; an observed zero is not proof of biological absence.",
    "Diet, weight change, and gastrointestinal symptoms after exposure may be mediators; covariate selection defines the estimand.")
  if (is.null(cfg$group_reference)) notes <- c(notes, "Single-cohort changes cannot separate exposure associations from secular time trends.")
  if (length(cfg$sensitivity_covariates)) notes <- c(notes, "Sensitivity covariates change the conditioning set; the sensitivity model is not a mediation or direct-effect analysis.")
  base_n <- tapply(as.character(md$phase) == "baseline", md$participant_id, sum)
  if (any(base_n < cfg$min_baseline_samples)) notes <- c(notes, "Some participants lack sufficient baseline samples; recovery is marked unavailable for them.")
  if (cfg$correlation == "independent") notes <- c(notes, "Residuals are conditionally independent given the random intercept; consider car1 for serial dependence.")
  batches <- intersect(c("assay", "body_site", "sequencing_batch"), names(md))
  for (nm in batches) if (length(unique(md[[nm]])) > 1)
    notes <- c(notes, paste("Multiple", nm, "values detected. Assay/body-site compatibility and batch adjustment require study-specific review."))
  bacteria <- read_abundance(cfg$bacteria, md, cfg, "bacteria")
  phages <- if (!is.null(cfg$phages)) read_abundance(cfg$phages, md, cfg, "phages") else NULL
  links <- NULL
  if (!is.null(cfg$host_links)) {
    links <- read_table(cfg$host_links)
    if (!all(c("phage_id","host_id","confidence","evidence") %in% names(links)))
      cb_stop("host_links requires phage_id, host_id, confidence, evidence.")
    if (anyNA(links[,c("phage_id","host_id","confidence","evidence")]) || anyDuplicated(links[,c("phage_id","host_id")]))
      cb_stop("Host links must be nonmissing and unique by phage/host pair; combine evidence explicitly.")
    if (any(!links$confidence %in% c("high","medium","low"))) cb_stop("Host confidence must be high, medium, or low.")
    if (any(!links$phage_id %in% colnames(phages$clr)) || any(!links$host_id %in% colnames(bacteria$clr)))
      cb_stop("Host links refer to missing or all-zero features. Reconcile identifiers explicitly.")
  }
  structure(list(config = cfg, metadata = md, bacteria = bacteria, phages = phages, links = links,
    notes = unique(notes)), class = "chronomicrobiome_study")
}
