chronomicrobiome_cli <- function(args = commandArgs(trailingOnly=TRUE)) {
  help <- paste(
    "Chrono\u00b5Biome 0.1.0 \u2014 longitudinal human microbiome and exposure analysis",
    "",
    "Usage:",
    "  chronomicrobiome doctor",
    "  chronomicrobiome demo --output PATH [--run]",
    "  chronomicrobiome validate --config study.yaml",
    "  chronomicrobiome run --config study.yaml [--output PATH]",
    "",
    "The demo is entirely synthetic. Output folders must be new or empty.",
    "R dependencies: nlme, ggplot2, yaml, jsonlite. Quarto is optional.",sep="\n")
  if (!length(args) || args[1] %in% c("help","--help","-h")) { cat(help,"\n"); return(invisible(NULL)) }
  if (identical(args[1],"--version")) { cat("chronomicrobiome 0.1.0\n"); return(invisible(NULL)) }
  command <- args[1]; flags <- list(); i <- 2L
  allowed <- switch(command,run=c("--config","--output"),validate="--config",demo=c("--output","--run"),doctor=character(),NULL)
  if (is.null(allowed)) cb_stop("Unknown command: ",command)
  while (i <= length(args)) {
    flag <- args[i]
    if (!flag %in% allowed || !is.null(flags[[flag]])) cb_stop("Unknown or repeated option: ",flag)
    if (flag=="--run") { flags[[flag]] <- TRUE; i <- i+1L } else {
      if (i==length(args) || startsWith(args[i+1L],"--")) cb_stop("Missing value for ",flag)
      flags[[flag]] <- args[i+1L]; i <- i+2L
    }
  }
  if (command=="doctor") {
    cat(R.version.string,"\n")
    for (p in c("nlme","ggplot2","yaml","jsonlite")) cat(p,": ",
      if (requireNamespace(p,quietly=TRUE)) as.character(utils::packageVersion(p)) else "MISSING", "\n",sep="")
    cat("Quarto: ",if (nzchar(find_quarto())) find_quarto() else "not found; built-in HTML available","\n",sep="")
    return(invisible(NULL))
  }
  if (command=="demo") {
    if (is.null(flags[["--output"]])) cb_stop("demo requires --output.")
    config <- make_demo(flags[["--output"]])
    cat("Synthetic study: ",config,"\n",sep="")
    if (isTRUE(flags[["--run"]])) return(invisible(run_study(config)))
    return(invisible(config))
  }
  if (is.null(flags[["--config"]])) cb_stop(command," requires --config.")
  if (command=="validate") {
    study <- validate_study(flags[["--config"]])
    cat("VALID: ",nrow(study$metadata)," samples; ",nlevels(study$metadata$participant_id)," participants\n",sep="")
    cat(paste0("- ",study$notes,collapse="\n"),"\n")
    return(invisible(study))
  }
  invisible(run_study(flags[["--config"]],flags[["--output"]]))
}
