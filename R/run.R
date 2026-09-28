write_tsv <- function(x,path) utils::write.table(x,path,sep="\t",quote=TRUE,row.names=FALSE,na="NA")

run_study <- function(config, output = NULL) {
  study <- validate_study(config)
  cfg <- study$config
  if (!is.null(output)) cfg$output <- path.expand(output)
  target <- cfg$output
  if (file.exists(target) && (!dir.exists(target) || length(list.files(target,all.files=TRUE,no..=TRUE))))
    cb_stop("Output must be new or empty; existing results are never overwritten: ", target)
  dir.create(target,recursive=TRUE,showWarnings=FALSE)
  target <- normalizePath(target,mustWork=TRUE)
  cfg$output <- target; study$config <- cfg
  dir.create(file.path(target,"tables")); dir.create(file.path(target,"figures"))
  started <- format(Sys.time(),tz="UTC",usetz=TRUE)
  writeLines(c("status: running",paste("started:",started)),file.path(target,"run-status.txt"))
  complete <- FALSE
  on.exit(if (!complete) writeLines("status: failed; partial outputs may exist",file.path(target,"run-status.txt")))
  message("Fitting longitudinal models...")
  fits <- fit_study(study)
  message("Computing persistence, return toward baseline, and host associations...")
  ecology <- ecology_study(study)
  hosts <- host_associations(study,fits)
  audit <- rbind(study$bacteria$audit,if (!is.null(study$phages)) study$phages$audit)
  tables <- c(list(coefficients=fits$coefficients,
    contrasts=fits$coefficients[fits$coefficients$primary_contrast,],model_diagnostics=fits$diagnostics),
    list(community_distances=ecology$distances,recovery=ecology$recovery,persistence=ecology$persistence,
      host_associations=hosts,feature_audit=audit))
  residual_rows <- list()
  for (key in names(fits$models)) for (feature in names(fits$models[[key]]$residuals)) {
    residual_rows[[paste(key,feature)]] <- data.frame(model=key,feature=feature,
      sample_id=study$metadata$sample_id,participant_id=study$metadata$participant_id,
      time=study$metadata$time,normalized_residual=fits$models[[key]]$residuals[[feature]],
      fitted_clr=fits$models[[key]]$fitted[[feature]])
  }
  tables$model_residuals <- do.call(rbind,residual_rows)
  for (nm in names(tables)) write_tsv(tables[[nm]],file.path(target,"tables",paste0(nm,".tsv")))
  write_tsv(study$metadata,file.path(target,"tables","metadata_used.tsv"))
  yaml::write_yaml(cfg[setdiff(names(cfg),"config_file")],file.path(target,"config-resolved.yaml"))
  utils::capture.output(utils::sessionInfo(),file=file.path(target,"session-info.txt"))
  message("Drawing figures and rendering the offline report...")
  figs <- plot_study(study,fits,ecology,target)
  report <- write_report(study,fits,ecology,hosts,figs,target)
  inputs <- unlist(cfg[c("metadata","bacteria","phages","host_links","config_file")],use.names=TRUE)
  packages <- c("nlme","ggplot2","yaml","jsonlite")
  provenance <- list(software="chronomicrobiome",version="0.1.0",status="complete",started_utc=started,
    finished_utc=format(Sys.time(),tz="UTC",usetz=TRUE),synthetic=cfg$synthetic,
    input_md5=as.list(setNames(as.character(tools::md5sum(inputs)),names(inputs))),input_paths=as.list(inputs),
    packages=as.list(setNames(vapply(packages,function(p)as.character(utils::packageVersion(p)),character(1)),packages)),
    report_engine=report$engine,seed=cfg$seed,notes=study$notes)
  jsonlite::write_json(provenance,file.path(target,"provenance.json"),pretty=TRUE,auto_unbox=TRUE)
  saveRDS(list(config=cfg,tables=tables,notes=study$notes),file.path(target,"analysis.rds"))
  writeLines(c("status: complete",paste("report_engine:",report$engine)),file.path(target,"run-status.txt"))
  complete <- TRUE
  message("Report: ",report$path)
  invisible(list(output=target,report=report$path,tables=tables,provenance=provenance))
}
