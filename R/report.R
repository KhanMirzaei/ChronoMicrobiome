utils::globalVariables(c("time","clr","participant_id","group","feature","phase","modality",
  "bray_curtis_to_baseline","estimate","ci_lower","ci_upper","label","residual","fitted","sample"))

html_escape <- function(x) {
  x <- as.character(x)
  for (pair in list(c("&","&amp;"),c("<","&lt;"),c(">","&gt;"),c('"',"&quot;"),c("'","&#39;")))
    x <- gsub(pair[1],pair[2],x,fixed=TRUE)
  x
}

html_table <- function(x, n = 20L) {
  if (is.null(x) || !nrow(x)) return('<p class="muted">No eligible results for this section.</p>')
  xx <- utils::head(x,n)
  cells <- lapply(xx,function(v) {
    if (is.numeric(v)) v <- ifelse(is.na(v), "Not available", format(signif(v,4),trim=TRUE))
    v[is.na(v)] <- "Not available"
    html_escape(v)
  })
  rows <- vapply(seq_len(nrow(xx)), function(i)
    paste0("<tr>",paste0("<td>",vapply(cells, `[`, character(1),i),"</td>",collapse=""),"</tr>"),character(1))
  paste0('<div class="table-wrap"><table><thead><tr>',paste0("<th>",html_escape(names(xx)),"</th>",collapse=""),
    "</tr></thead><tbody>",paste(rows,collapse=""),"</tbody></table></div>",
    if (nrow(x)>n) paste0('<p class="muted">Showing ',n," of ",nrow(x)," rows. Complete results are in the TSV downloads.</p>") else "")
}

plot_study <- function(study, fits, ecology, output) {
  cfg <- study$config; md <- study$metadata
  cap <- if (cfg$synthetic) "SYNTHETIC DEMONSTRATION \u2014 not patient data or evidence about semaglutide" else "Exploratory associations; review model assumptions and study design."
  theme <- ggplot2::theme_minimal(base_size=12) + ggplot2::theme(
    panel.grid.minor=ggplot2::element_blank(), plot.title=ggplot2::element_text(face="bold",colour="#15333b"),
    strip.text=ggplot2::element_text(face="bold"), legend.position="bottom",
    plot.caption=ggplot2::element_text(colour="#7c4b35",size=9,hjust=0), plot.background=ggplot2::element_rect(fill="white",colour=NA))
  save <- function(p, name, height=6) {
    path <- file.path(output,"figures",paste0(name,".png"))
    ggplot2::ggsave(path, p + theme, width=11,height=height,dpi=150,bg="white")
    path
  }
  figs <- list()
  p <- ggplot2::ggplot(md,ggplot2::aes(x=time,y=participant_id,colour=phase,group=participant_id)) +
    ggplot2::geom_line(colour="#d5deda") + ggplot2::geom_point(size=1.9) +
    ggplot2::facet_grid(group ~ .,scales="free_y",space="free_y") +
    ggplot2::labs(title="Sampling and exposure windows",x="Days relative to study index event",y="Participant",caption=cap)
  figs$timeline <- save(p,"01_sampling",7)
  trajectories <- list()
  for (mod in c("bacteria", if (!is.null(study$phages)) "phages")) {
    for (f in utils::head(study[[mod]]$keep,3)) trajectories[[paste(mod,f)]] <-
      data.frame(time=md$time,participant_id=md$participant_id,group=md$group,
        feature=paste(mod,f,sep=" \u00b7 "),clr=study[[mod]]$clr[,f])
  }
  td <- do.call(rbind,trajectories)
  p <- ggplot2::ggplot(td,ggplot2::aes(x=time,y=clr,group=participant_id,colour=group)) +
    ggplot2::geom_line(alpha=0.25,linewidth=0.45) + ggplot2::geom_point(alpha=0.35,size=0.7) +
    ggplot2::facet_wrap(~ feature,scales="free_y",ncol=2) +
    ggplot2::scale_colour_manual(values=c("#1d807c","#c46b42")) +
    ggplot2::labs(title="Individual microbial trajectories",subtitle="First three eligible features per modality; natural-log CLR scale",
      x="Days relative to study index event",y="Centered log-ratio abundance",caption=cap)
  figs$trajectories <- save(p,"02_trajectories",8)
  contrasts <- fits$coefficients[fits$coefficients$primary_contrast & fits$coefficients$specification=="primary",]
  if (nrow(contrasts)) {
    fd <- utils::head(contrasts[order(-abs(contrasts$estimate)),],12)
    short_term <- sub("^.*phase", "", fd$term)
    fd$label <- factor(paste(fd$feature,short_term,sep=" / "), levels=rev(unique(paste(fd$feature,short_term,sep=" / "))))
    p <- ggplot2::ggplot(fd,ggplot2::aes(x=estimate,y=label,colour=modality)) +
      ggplot2::geom_vline(xintercept=0,linetype=2,colour="#a9b6b2") +
      ggplot2::geom_segment(ggplot2::aes(x=ci_lower,xend=ci_upper,yend=label),linewidth=0.7) +
      ggplot2::geom_point(size=2.5) + ggplot2::scale_colour_manual(values=c("#1d807c","#c46b42")) +
      ggplot2::labs(title="Exposure-associated contrasts",subtitle="Up to 12 largest absolute estimates; unadjusted 95% confidence intervals",
        x=if (is.null(cfg$group_reference)) "Change from baseline on the CLR scale" else "Difference in change on the CLR scale",y=NULL,caption=cap)
    figs$effects <- save(p,"03_effects",7)
  }
  dd <- ecology$distances[is.finite(ecology$distances$bray_curtis_to_baseline),]
  if (nrow(dd)) {
    p <- ggplot2::ggplot(dd,ggplot2::aes(x=time,y=bray_curtis_to_baseline,group=participant_id,colour=group)) +
      ggplot2::geom_line(alpha=0.4) + ggplot2::geom_point(alpha=0.5,size=1) +
      ggplot2::geom_hline(yintercept=cfg$recovery_threshold,linetype=2,colour="#993d2b") +
      ggplot2::facet_wrap(~modality,ncol=1) + ggplot2::scale_colour_manual(values=c("#1d807c","#c46b42")) +
      ggplot2::labs(title="Community distance from personal baseline",subtitle="Dashed line: user-defined descriptive return threshold",
        x="Days relative to study index event",y="Bray\u2013Curtis distance",caption=cap)
    figs$recovery <- save(p,"04_recovery",7)
  }
  rd <- list()
  for (mod in c("bacteria",if (!is.null(study$phages)) "phages")) for (f in utils::head(names(fits$models[[mod]]$residuals),2))
    rd[[paste(mod,f)]] <- data.frame(feature=paste(mod,f),residual=fits$models[[mod]]$residuals[[f]],
      fitted=fits$models[[mod]]$fitted[[f]])
  if (length(rd)) {
    rd <- do.call(rbind,rd)
    p <- ggplot2::ggplot(rd,ggplot2::aes(x=fitted,y=residual)) + ggplot2::geom_point(alpha=0.35,colour="#1d807c") +
      ggplot2::geom_hline(yintercept=0,linetype=2) + ggplot2::facet_wrap(~feature,scales="free_x") +
      ggplot2::labs(title="Model diagnostics: residuals versus fitted",subtitle="First two successfully fitted features per modality",
        x="Fitted CLR",y="Normalized residual",caption=cap)
    figs$diagnostics <- save(p,"05_residuals")
    p <- ggplot2::ggplot(rd,ggplot2::aes(sample=residual)) + ggplot2::stat_qq(colour="#1d807c",alpha=0.5) +
      ggplot2::stat_qq_line(colour="#c46b42") + ggplot2::facet_wrap(~feature) +
      ggplot2::labs(title="Model diagnostics: normal Q\u2013Q plots",x="Theoretical quantile",y="Normalized residual quantile",caption=cap)
    figs$qq <- save(p,"06_qq")
  }
  figs
}

find_quarto <- function() {
  candidates <- c(Sys.getenv("CHRONOMICROBIOME_QUARTO"), Sys.which("quarto"),
    "/Applications/RStudio.app/Contents/Resources/app/quarto/bin/quarto",
    "/Applications/RStudio.app/Contents/Resources/quarto/bin/quarto")
  candidates <- candidates[nzchar(candidates) & file.exists(candidates)]
  if (length(candidates)) candidates[1] else ""
}

write_report <- function(study, fits, ecology, hosts, figures, output) {
  cfg <- study$config
  image_tag <- function(key, alt) {
    path <- figures[[key]]
    if (is.null(path)) return('<p class="muted">No eligible data for this plot.</p>')
    bytes <- readBin(path,"raw",n=file.info(path)$size)
    paste0('<img class="figure" alt="',html_escape(alt),'" src="data:image/png;base64,',jsonlite::base64_enc(bytes),'">')
  }
  contrasts <- fits$coefficients[fits$coefficients$primary_contrast,]
  primary <- contrasts[contrasts$specification=="primary",]
  primary <- primary[order(primary$q_value),]
  sensitivity <- contrasts[contrasts$specification=="sensitivity",]
  display_cols <- c("modality","feature","term","estimate","ci_lower","ci_upper","q_value")
  recovery_counts <- stats::aggregate(list(participants=ecology$recovery$participant_id),
    ecology$recovery[,c("modality","group","status")],length)
  eligible <- hosts[hosts$included_primary,,drop=FALSE]
  high <- hosts[hosts$included_high_only,,drop=FALSE]
  formula <- paste(deparse(model_formula(cfg)),collapse=" ")
  css <- paste(readLines(cb_file("report","report.css"),warn=FALSE),collapse="\n")
  banner <- if (cfg$synthetic) '<div class="notice"><strong>Synthetic demonstration.</strong> Every participant, abundance, exposure pattern, and result on this page is simulated. These results do not describe the effects of Ozempic or semaglutide in humans.</div>' else
    '<div class="notice">Exploratory cohort analysis. Estimates describe associations under the specified model; they do not establish causal medication effects.</div>'
  estimand <- if (is.null(cfg$group_reference)) "Phase coefficients estimate within-cohort changes relative to baseline. Without a comparator, exposure and secular time effects cannot be separated." else
    paste0("Interaction coefficients estimate the difference in baseline-to-phase change between ",html_escape(levels(study$metadata$group)[2]),
      " and ",html_escape(cfg$group_reference),". Causal interpretation requires a suitable design and additional assumptions; it is not supplied by this model.")
  body <- paste0('<main><header><p class="eyebrow">CHRONO\u00b5BIOME / RESEARCH PREVIEW 0.1.0</p><h1>',html_escape(cfg$title),
    '</h1><p class="subtitle">Human microbiome dynamics and responses to exposures</p></header>',banner,
    '<nav><a href="#design">Study design</a><a href="#effects">Exposure associations</a><a href="#recovery">Recovery</a><a href="#hosts">Phage\u2013host links</a><a href="#methods">Methods & quality</a></nav>',
    '<div class="metrics"><div><b>',nlevels(study$metadata$participant_id),'</b><span>participants</span></div><div><b>',nrow(study$metadata),
    '</b><span>samples</span></div><div><b>',length(study$bacteria$keep),'</b><span>bacterial features modeled</span></div><div><b>',
    if (is.null(study$phages)) 0 else length(study$phages$keep),'</b><span>phage features modeled</span></div></div>',
    '<section id="design"><h2>Study design and individual trajectories</h2><p>',estimand,'</p>',
    image_tag("timeline","Participant sampling timelines"),image_tag("trajectories","Individual CLR abundance trajectories"),'</section>',
    '<section id="effects"><h2>Exposure-associated changes</h2><p>Estimates use natural-log centered log-ratio (CLR) abundance. They are relative to the geometric mean of the supplied nonzero features within each modality, not absolute abundance changes or ordinary fold changes.</p>',
    image_tag("effects","Mixed-model exposure contrasts"),html_table(primary[,display_cols,drop=FALSE]),
    '<p class="muted">BH-adjusted p-values cover all planned feature \u00d7 nonbaseline-phase contrasts across both modalities, separately for each model specification. Confidence intervals are pointwise, not multiplicity-adjusted. Other coefficients are available in the download and are not assigned q-values.</p>',
    if (length(cfg$sensitivity_covariates)) paste0('<details><summary>Sensitivity model: additional conditioning variables</summary><p>Added variables: ',
      html_escape(paste(cfg$sensitivity_covariates,collapse=", ")),'. These may be downstream of treatment. This is a sensitivity analysis, not an estimate of a direct or mediated drug effect.</p>',
      html_table(sensitivity[,display_cols,drop=FALSE]),'</details>') else "",'</section>',
    '<section id="recovery"><h2>Persistence and return toward baseline</h2><p>Community distances use the complete supplied composition, before prevalence filtering. Baseline points use a leave-one-out reference. A descriptive return requires ',
    cfg$recovery_consecutive,' consecutive observed follow-up samples with distance \u2264 ',cfg$recovery_threshold,
    '. This configurable threshold is not a validated biological recovery criterion. Return is evaluated only after an observed departure during an exposure phase.</p>',
    image_tag("recovery","Bray-Curtis distance from personal baseline"),html_table(recovery_counts),
    '<p>First observed return and confirmation dates are sampling dates, not exact biological recovery times. Missing follow-up is not counted as failure to recover. Feature persistence is the fraction of observed samples above the configured detection threshold; it does not establish continuous colonization.</p></section>',
    '<section id="hosts"><h2>Phage\u2013host trajectory associations</h2><p>',nrow(eligible),' supplied links meet the primary confidence filter; ',nrow(high),
    ' meet the high-confidence-only filter. Confidence categories are supplied evidence labels, not calibrated probabilities.</p>',
    html_table(eligible[,setdiff(names(eligible),c("included_primary","included_high_only")),drop=FALSE]),
    '<details><summary>High-confidence-only sensitivity</summary>',html_table(high[,setdiff(names(high),c("included_primary","included_high_only")),drop=FALSE]),'</details>',
    '<p>These are means of within-participant Spearman correlations of model residuals. Percentile intervals resample participants conditional on fitted models; they do not include model-fitting or host-assignment uncertainty and are not adjusted for multiple comparisons. At least four samples per participant and six eligible participants per pair are required. Correlation does not establish infection, killing, or direction of causation. Lagged inference is not implemented in this release.</p></section>',
    '<section id="methods"><h2>Methods, diagnostics, and provenance</h2><pre>',html_escape(formula),
    '\nRandom intercept: participant_id\nResidual correlation: ',html_escape(cfg$correlation),
    '\nCLR pseudocount (relative-abundance units): ',cfg$pseudocount,'\nPrevalence threshold: ',cfg$min_prevalence,
    '\nDetection threshold (relative-abundance units): ',cfg$detection_threshold,'</pre>',
    '<ul>',paste0('<li>',html_escape(study$notes),'</li>',collapse=""),'</ul>',
    '<details open><summary>Model fitting status</summary>',html_table(fits$diagnostics),'</details>',
    image_tag("diagnostics","Residual versus fitted diagnostic plots"),image_tag("qq","Normal quantile diagnostic plots"),
    '<p>The exploratory CLR model is not a substitute for study-specific differential-abundance validation. Rare taxa, zero replacement, unmeasured confounding, dropout, and serial dependence can affect results. Review diagnostics for every modeled feature using the saved residual table, not just the example plots.</p>',
    '<h3>Complete outputs</h3><ul>',paste0('<li><a href="tables/',
      c("contrasts.tsv","coefficients.tsv","model_diagnostics.tsv","model_residuals.tsv","community_distances.tsv","recovery.tsv","persistence.tsv","host_associations.tsv","feature_audit.tsv"),
      '">',c("Exposure contrasts","All model coefficients","Model status","Model residuals","Community distances","Return-to-baseline status","Observed persistence","Host evidence associations","Feature filtering audit"),'</a></li>',collapse=""),
    '</ul><p><a href="provenance.json">Run provenance and input checksums</a> \u00b7 <a href="session-info.txt">R session information</a> \u00b7 <a href="config-resolved.yaml">Resolved configuration</a></p>',
    '<p class="muted">Chrono\u00b5Biome (chronomicrobiome) is not affiliated with the existing chronobiome.org research project. This release has synthetic-data checks but has not been externally validated on a real cohort.</p></section>',
    '<footer>Chrono\u00b5Biome 0.1.0 \u00b7 Reproducible, local, exploratory analysis</footer></main>')
  path <- file.path(output,"report.html")
  fallback <- paste0('<!doctype html><html lang="en"><head><meta charset="utf-8"><meta name="viewport" content="width=device-width,initial-scale=1"><title>',
    html_escape(cfg$title),'</title><style>',css,'</style></head><body>',body,'</body></html>')
  writeLines(fallback,path,useBytes=TRUE)
  engine <- "builtin"
  quarto <- find_quarto()
  if (cfg$report_engine == "quarto" && !nzchar(quarto)) cb_stop("Quarto requested but unavailable. Built-in report saved; install Quarto or choose auto/builtin.")
  if (cfg$report_engine != "builtin" && nzchar(quarto)) {
    file.copy(cb_file("report","report.css"),file.path(output,"report.css"),overwrite=TRUE)
    # Keep data-derived HTML out of executable Quarto/Markdown source.
    writeLines(body,file.path(output,"report-body.html"),useBytes=TRUE)
    file.copy(cb_file("report","template.qmd"),file.path(output,"report.qmd"),overwrite=TRUE)
    status <- tryCatch(system2(quarto,c("render",shQuote(file.path(output,"report.qmd")),"--to","html","--output","report.html"),
      stdout=file.path(output,"quarto.log"),stderr=file.path(output,"quarto-errors.log")),error=function(e) 1L)
    if (status == 0L && file.exists(path)) engine <- "quarto" else {
      writeLines(fallback,path,useBytes=TRUE)
      if (cfg$report_engine == "quarto") cb_stop("Quarto failed; see quarto-errors.log. A built-in report was saved.")
      warning("Quarto failed; using the built-in offline report. See quarto-errors.log.",call.=FALSE)
    }
  }
  list(path=path,engine=engine)
}
