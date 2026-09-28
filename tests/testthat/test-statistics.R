test_that("injected directions are recovered and BH covers the full family", {
  f <- fixture(); on.exit(unlink(f$root,recursive=TRUE))
  study <- validate_study(f$path)
  fits <- chronomicrobiome:::fit_study(study)
  x <- fits$coefficients
  select <- x$specification=="primary" & x$primary_contrast & grepl("phasetreatment",x$term)
  expect_gt(x$estimate[select & x$feature=="Bacterium_01"],0.6)
  expect_lt(x$estimate[select & x$feature=="Bacterium_02"],-0.5)
  idx <- x$specification=="primary" & x$primary_contrast
  planned <- sum(fits$diagnostics$specification=="primary") * 2
  expect_equal(x$q_value[idx],p.adjust(x$p_value[idx],"BH",n=planned))
  expect_true(all(is.na(x$q_value[!x$primary_contrast])))
  expect_true(all(x$ci_lower <= x$estimate & x$ci_upper >= x$estimate))
  h <- chronomicrobiome:::host_associations(study,fits)
  expect_equal(nrow(h),4L)
  expect_equal(sum(h$included_high_only),2L)
  expect_equal(sum(h$included_primary),3L)
  expect_true(all(h$n_participants==32L))
  expect_true(all(h$ci_lower >= -1 & h$ci_upper <= 1))
})

test_that("continuous-time correlation fits and flat responses fail transparently", {
  f <- fixture(); on.exit(unlink(f$root,recursive=TRUE))
  study <- validate_study(f$path)
  cfg <- study$config; cfg$correlation <- "car1"
  fit <- chronomicrobiome:::fit_feature(study$bacteria$clr[,1],study$metadata,cfg)
  expect_null(fit$error)
  expect_true(is.finite(fit$rho) && fit$rho >= 0 && fit$rho <= 1)
  flat <- chronomicrobiome:::fit_feature(rep(0,nrow(study$metadata)),study$metadata,cfg)
  expect_match(flat$error,"No response variation")
})

test_that("recovery requires observed departure and consecutive follow-up samples", {
  md <- data.frame(sample_id=paste0("S",1:8),participant_id=factor(rep("P1",8)),group="A",
    time=c(-14,-7,0,7,14,21,28,35),phase=c("baseline","baseline","treatment","treatment",rep("followup",4)))
  x <- c(0.48,0.52,0.9,0.85,0.52,0.51,0.53,0.55)
  data <- list(relative=cbind(x,1-x))
  cfg <- list(min_baseline_samples=2,recovery_threshold=0.1,recovery_consecutive=2)
  result <- chronomicrobiome:::community_recovery(data,md,cfg,"bacteria")
  expect_equal(result$distances$bray_curtis_to_baseline[1],0.04)
  expect_equal(result$recovery$status,"observed_return_to_threshold")
  expect_equal(result$recovery$first_observed_return_day,14)
  expect_equal(result$recovery$confirmation_day,21)
  data$relative[3:4,] <- 0.5
  expect_equal(chronomicrobiome:::community_recovery(data,md,cfg,"bacteria")$recovery$status,"no_observed_departure")
  cfg$min_baseline_samples <- 3
  expect_equal(chronomicrobiome:::community_recovery(data,md,cfg,"bacteria")$recovery$status,"insufficient_baseline")
})

test_that("seeded functions preserve RNG state", {
  set.seed(7); before <- .Random.seed
  f <- fixture(); on.exit(unlink(f$root,recursive=TRUE))
  expect_identical(.Random.seed,before)
})
