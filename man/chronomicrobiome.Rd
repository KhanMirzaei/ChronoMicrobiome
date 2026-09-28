\name{chronomicrobiome}
\alias{run_study}
\alias{validate_study}
\alias{make_demo}
\alias{chronomicrobiome_cli}
\title{Longitudinal Microbiome and Exposure Analysis}
\description{Validate abundance tables and repeated-measures metadata, fit exploratory CLR mixed models, describe return toward baseline, and report phage-host associations.}
\usage{
run_study(config, output = NULL)
validate_study(config)
make_demo(output, seed = 20260928L)
chronomicrobiome_cli(args = commandArgs(trailingOnly = TRUE))
}
\arguments{
\item{config}{Path to a YAML study configuration. Input paths are resolved relative to this file.}
\item{output}{A new or empty output directory. For \code{run_study}, overrides the configured path and is relative to the working directory.}
\item{seed}{Integer seed used to generate synthetic data without changing the caller's random-number state.}
\item{args}{Character vector of command-line arguments.}
}
\details{Metadata requires sample_id, participant_id, numeric time in days, and phase. Optional two-group designs require a constant group per participant and group_reference in YAML. Wide abundance tables have sample_id first and at least two feature columns. All modalities must exactly match metadata sample IDs. Counts or proportions are supported, specified by abundance_scale. The first phase must be baseline. Optional followup enables descriptive recovery status.

Phase coefficients in a single cohort describe changes from baseline. Group-by-phase interaction coefficients describe differences in changes between groups. These are exploratory relative-abundance associations, not causal drug effects. Primary multiple-testing correction covers planned feature-by-phase contrasts across modalities. Confidence intervals are pointwise.

Missing required metadata and incompatible identifiers are rejected explicitly. Reports include diagnostics and input hashes. This version supports one ordered exposure episode, with no lagged, mediation, dose-response, or causal inference. See the source README and methods documentation for all restrictions.}
\value{\code{validate_study} returns a validated study list. \code{run_study} invisibly returns paths, tables, and provenance. \code{make_demo} invisibly returns the synthetic YAML path. \code{chronomicrobiome_cli} returns the command result invisibly.}
\examples{
directory <- tempfile("chronomicrobiome-example-")
config <- make_demo(directory)
study <- validate_study(config)
nrow(study$metadata)
unlink(directory, recursive = TRUE)
}
