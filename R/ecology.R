community_recovery <- function(data, md, cfg, label) {
  distances <- list(); returns <- list()
  for (id in levels(md$participant_id)) {
    idx <- which(md$participant_id == id)
    base <- idx[as.character(md$phase[idx]) == "baseline"]
    enough <- length(base) >= cfg$min_baseline_samples
    vals <- rep(NA_real_, length(idx))
    if (enough) for (j in seq_along(idx)) {
      ref <- setdiff(base, idx[j])
      if (length(ref)) vals[j] <- sum(abs(data$relative[idx[j],] - colMeans(data$relative[ref,,drop=FALSE]))) / 2
    }
    distances[[id]] <- data.frame(modality=label, sample_id=md$sample_id[idx], participant_id=id,
      group=as.character(md$group[idx]), time=md$time[idx], phase=as.character(md$phase[idx]),
      bray_curtis_to_baseline=vals, n_baseline=length(base))
    follow <- which(as.character(md$phase[idx]) == "followup")
    exposed <- which(!as.character(md$phase[idx]) %in% c("baseline","followup"))
    status <- if (!enough) "insufficient_baseline" else if (!length(exposed)) "insufficient_exposure_samples" else
      if (!length(follow)) "no_followup" else if (all(vals[exposed] <= cfg$recovery_threshold)) "no_observed_departure" else "not_observed_by_last_sample"
    first <- confirmation <- lower <- NA_real_
    if (status == "not_observed_by_last_sample") {
      within <- vals[follow] <= cfg$recovery_threshold
      k <- cfg$recovery_consecutive
      if (length(follow) >= k) for (j in seq_len(length(follow)-k+1L)) {
        if (all(within[j:(j+k-1L)])) {
          status <- "observed_return_to_threshold"
          first <- md$time[idx[follow[j]]]
          confirmation <- md$time[idx[follow[j+k-1L]]]
          lower <- if (j > 1) md$time[idx[follow[j-1L]]] else NA_real_
          break
        }
      }
    }
    returns[[id]] <- data.frame(modality=label, participant_id=id, group=as.character(md$group[idx[1]]),
      status=status, n_baseline=length(base), n_followup=length(follow),
      first_observed_return_day=first, confirmation_day=confirmation,
      previous_followup_day=lower, last_observed_day=max(md$time[idx]), threshold=cfg$recovery_threshold,
      consecutive_samples=cfg$recovery_consecutive)
  }
  list(distances=do.call(rbind,distances), recovery=do.call(rbind,returns))
}

feature_persistence <- function(data, md, cfg, label) {
  rows <- lapply(levels(md$participant_id), function(id) {
    idx <- which(md$participant_id == id)
    detected <- data$relative[idx,,drop=FALSE] > cfg$detection_threshold
    data.frame(modality=label, participant_id=id, feature=colnames(detected),
      n_samples=length(idx), n_detected=colSums(detected), observed_detection_fraction=colMeans(detected),
      first_detected_day=vapply(seq_len(ncol(detected)), function(j) if (any(detected[,j])) min(md$time[idx][detected[,j]]) else NA_real_, numeric(1)),
      last_detected_day=vapply(seq_len(ncol(detected)), function(j) if (any(detected[,j])) max(md$time[idx][detected[,j]]) else NA_real_, numeric(1)))
  })
  do.call(rbind,rows)
}

ecology_study <- function(study) {
  ec <- list(); pers <- list()
  for (mod in c("bacteria", if (!is.null(study$phages)) "phages")) {
    ec[[mod]] <- community_recovery(study[[mod]], study$metadata, study$config, mod)
    pers[[mod]] <- feature_persistence(study[[mod]], study$metadata, study$config, mod)
  }
  list(distances=do.call(rbind,lapply(ec, `[[`, "distances")),
    recovery=do.call(rbind,lapply(ec, `[[`, "recovery")), persistence=do.call(rbind,pers))
}
