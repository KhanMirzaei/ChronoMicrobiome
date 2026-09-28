make_demo <- function(output, seed = 20260928L) {
  if (file.exists(output) && (!dir.exists(output) || length(list.files(output,all.files=TRUE,no..=TRUE))))
    cb_stop("Demo destination must be new or empty: ",output)
  dir.create(output,recursive=TRUE,showWarnings=FALSE)
  output <- normalizePath(output,mustWork=TRUE)
  with_cb_seed(seed,function() {
    days <- c(-14,-7,0,14,28,56,84,112,140,168)
    n <- 32L
    rows <- list(); b_rows <- list(); p_rows <- list()
    bn <- sprintf("Bacterium_%02d",seq_len(12)); pn <- sprintf("Phage_%02d",seq_len(8))
    for (i in seq_len(n)) {
      treated <- i > n/2
      group <- if (treated) "semaglutide" else "comparator"
      age <- round(stats::runif(1,30,65)); bmi <- round(stats::runif(1,28,40),1)
      metformin <- if (i %% 3 == 0) "yes" else "no"
      b_person <- stats::rnorm(12,0,0.45); p_person <- stats::rnorm(8,0,0.55)
      for (j in seq_along(days)) {
        time <- days[j]
        phase <- if (time < 0) "baseline" else if (time <= 84) "treatment" else "followup"
        active <- as.numeric(treated && phase == "treatment")
        response_strength <- if (time < 0) 0 else if (time <= 84) 1 else exp(-(time-84)/35)
        effect <- as.numeric(treated) * response_strength
        weight <- if (time<0) 0 else -as.numeric(treated)*5*(1-exp(-(time+1)/35)) + stats::rnorm(1,0,0.7)
        fiber <- if (time<0) 0 else as.numeric(treated)*2 + stats::rnorm(1,0,1.5)
        shared <- stats::rnorm(1,0,0.25)
        blog <- b_person + stats::rnorm(12,0,0.2)
        blog[1] <- blog[1] + 1.1*effect + shared
        blog[2] <- blog[2] - 0.85*effect
        blog[3] <- blog[3] + (age-45)*0.006
        plog <- p_person + stats::rnorm(8,0,0.25)
        plog[1] <- plog[1] + 0.8*effect + shared
        plog[2] <- plog[2] - 0.6*effect
        bp <- exp(blog-max(blog)); bp <- bp/sum(bp)
        pp <- exp(plog-max(plog)); pp <- pp/sum(pp)
        if (stats::runif(1)<0.65) { bp[12] <- 0; bp <- bp/sum(bp) }
        if (stats::runif(1)<0.6) { pp[8] <- 0; pp <- pp/sum(pp) }
        sid <- sprintf("SYN_P%02d_T%02d",i,j)
        rows[[sid]] <- data.frame(sample_id=sid,participant_id=sprintf("SYN_P%02d",i),time=time,
          phase=phase,group=group,age_baseline=age,bmi_baseline=bmi,metformin_baseline=metformin,
          exposure_on=active,weight_change_pct=round(weight,3),fiber_change_g_day=round(fiber,3),
          body_site="stool",assay="bulk_shotgun",synthetic=TRUE)
        b_rows[[sid]] <- stats::setNames(as.numeric(stats::rmultinom(1,50000,bp)),bn)
        p_rows[[sid]] <- stats::setNames(as.numeric(stats::rmultinom(1,8000,pp)),pn)
      }
    }
    md <- do.call(rbind,rows); rownames(md) <- NULL
    b <- data.frame(sample_id=names(b_rows),do.call(rbind,b_rows),check.names=FALSE,row.names=NULL)
    p <- data.frame(sample_id=names(p_rows),do.call(rbind,p_rows),check.names=FALSE,row.names=NULL)
    links <- data.frame(phage_id=pn[1:4],host_id=bn[1:4],confidence=c("high","high","medium","low"),
      evidence=c("Synthetic integration evidence","Synthetic CRISPR evidence","Synthetic genus-level prediction","Synthetic weak prediction"))
    write_tsv(md,file.path(output,"metadata.tsv")); write_tsv(b,file.path(output,"bacteria.tsv"))
    write_tsv(p,file.path(output,"phages.tsv")); write_tsv(links,file.path(output,"host_links.tsv"))
    cfg <- list(title="Semaglutide exposure: synthetic study",synthetic=TRUE,
      metadata="metadata.tsv",bacteria="bacteria.tsv",phages="phages.tsv",host_links="host_links.tsv",
      output="results",abundance_scale="counts",phase_levels=c("baseline","treatment","followup"),
      group_reference="comparator",covariates=c("age_baseline","bmi_baseline","metformin_baseline"),
      sensitivity_covariates=c("weight_change_pct","fiber_change_g_day"),
      min_prevalence=0.1,pseudocount=1e-6,correlation="independent",time_unit_days=7,
      min_baseline_samples=2L,recovery_threshold=0.15,recovery_consecutive=2L,detection_threshold=0,
      host_confidence=c("high","medium"),bootstrap_replicates=500L,seed=as.integer(seed),report_engine="auto")
    yaml::write_yaml(cfg,file.path(output,"study.yaml"))
    writeLines(c("FULLY SYNTHETIC EXAMPLE \u2014 NO PATIENT DATA",
      "32 simulated participants, 10 observations each; comparator and semaglutide labels are illustrative.",
      "All feature identities and abundance changes are invented; they are not findings about any drug.",
      "The illustrative exposure and follow-up schedule is not prescribing or discontinuation guidance.",
      "Known simulation signals: Bacterium_01 increases; Bacterium_02 decreases in the exposure group.",
      "Phage_01 shares a simulated residual driver with Bacterium_01. Feature 12/8 include synthetic zeros.",
      "Phage abundance is separately normalized. Cross-modality absolute phage/host ratios are not inferred."),
      file.path(output,"README.txt"))
  })
  invisible(file.path(output,"study.yaml"))
}
