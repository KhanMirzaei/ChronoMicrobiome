# ChronoµBiome

**Human microbiome dynamics and responses to exposures.**

ChronoµBiome 0.1.0 is a working R package and command-line workflow for processed
longitudinal bacterial and optional phage abundance tables. It fits exploratory
exposure-associated models, describes observed persistence and return toward
baseline, and summarizes supplied phage–host links against repeated measurements.

**Research preview:** tested on synthetic data; the included application is entirely simulated and contains no evidence about actual medication effects.

The display name is **ChronoµBiome** (pronounced “Chrono-micro-biome”). The R
package and command use the ASCII name `chronomicrobiome`. This project is
independent of the unrelated [Chronobiome research project](https://chronobiome.org/).

[![R package checks](https://github.com/KhanMirzaei/ChronoMicrobiome/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/KhanMirzaei/ChronoMicrobiome/actions/workflows/R-CMD-check.yaml)

**[Explore the synthetic demo online](https://khanmirzaei.github.io/ChronoMicrobiome/)**

## Get the source

```bash
git clone https://github.com/KhanMirzaei/ChronoMicrobiome.git
cd ChronoMicrobiome
```

Or install the R package directly from GitHub:

```r
install.packages("remotes")
remotes::install_github("KhanMirzaei/ChronoMicrobiome")
```

## Quick start

Requirements: R 4.1 or later and four CRAN packages. Quarto is optional.

```r
install.packages(c("nlme", "ggplot2", "yaml", "jsonlite"))
```

From the extracted ChronoµBiome source directory:

```bash
Rscript bin/chronomicrobiome doctor
Rscript bin/chronomicrobiome demo --output my-synthetic-study --run
```

Open `my-synthetic-study/results/report.html`. The report includes embedded
figures and needs no internet connection. Keep the results folder together to
preserve links to downloadable tables and provenance.

The source launcher works without installing ChronoµBiome itself. On macOS or Linux,
`./bin/chronomicrobiome` is equivalent to `Rscript bin/chronomicrobiome`. Commands and file
paths containing spaces are supported.

An existing synthetic study is also included:

```bash
Rscript bin/chronomicrobiome validate --config inst/extdata/semaglutide/study.yaml
Rscript bin/chronomicrobiome run --config inst/extdata/semaglutide/study.yaml --output my-results
```

Output directories must be new or empty. Existing results are never overwritten.
Config-relative input/output paths are resolved against the YAML file. A CLI
`--output` override is relative to your current working directory.

## Install as an R package

```bash
R CMD INSTALL .
```

```r
library(chronomicrobiome)
validate_study("study.yaml")
result <- run_study("study.yaml", output = "results")
result$report
```

The source launcher is in `bin/`; installing the R package does not automatically
put a shell command on PATH. Within R, `chronomicrobiome_cli(c("doctor"))` is also available.

## What version 0.1 does

1. Validates metadata, identifiers, counts/proportions, design rank, and host links.
2. Normalizes each modality separately and computes natural-log CLR abundances.
3. Fits random-intercept linear mixed models using `nlme`, with optional continuous-
   time AR(1) residual correlation for irregularly spaced visits.
4. Reports group-by-phase contrasts (or phase changes in a single cohort),
   pointwise confidence intervals, and BH-adjusted p-values.
5. Summarizes observed feature persistence and community distances from personal
   baseline; records observed return, unavailable follow-up, and inadequate baseline.
6. Summarizes within-person residual correlations for supplied phage–host pairs,
   with a participant bootstrap and confidence-filter sensitivity.
7. Saves TSV tables, figures, resolved configuration, MD5 input fingerprints,
   R session/package versions, an RDS results object, and an HTML report.

This release supports **one ordered exposure episode** and either one cohort or
two fixed study groups. Medication, antibiotic, and diet studies use the same
configuration interface. In v0.1, dose and adherence can be recorded in metadata
but do not receive a dedicated dose-response analysis.

## Input format

Use tab-separated `.tsv` or comma-separated `.csv` files.

**Metadata: one row per biological sample.**

| Column | Required | Meaning |
|---|---|---|
| `sample_id` | Yes | Unique, nonmissing sample identifier |
| `participant_id` | Yes | Person sampled repeatedly |
| `time` | Yes | Numeric days relative to the index event; same definition across people |
| `phase` | Yes | Ordered window such as baseline, treatment, followup |
| `group` | For comparisons | Exactly two groups, constant within each person |
| Covariate columns | If configured | Numeric values or categorical labels, fully observed |

Control participants need an explicitly defined corresponding index date and
equivalent windows. ChronoµBiome does not infer exposure dates or phases from drug
names. For a recovery analysis, label post-exposure observation windows `followup`;
do not use this label to imply a medication stopped if treatment continued.

**Abundance: one row per sample, features as columns.**

```text
sample_id    Bacterium_01    Bacterium_02    Bacterium_03
S01          120             300             80
S02          150             250             100
```

The first column must be `sample_id`. Each modality needs at least two nonzero
features and must match the metadata sample set exactly. Missing observations
must be resolved before analysis; the tool never substitutes zeros for a missing
sample. Technical replicates must be resolved explicitly before input.

Set `abundance_scale: counts` for integer counts or `relative` for proportions
summing to 1 (not percentages). The scale applies to both modalities. Supply a
consistent, sufficiently broad feature table rather than a hand-picked pair of
taxa: CLR interpretation depends on the supplied reference composition. All-zero
features are excluded; other filtering is audited.

**Optional host links:** `phage_id`, `host_id`, `confidence`, `evidence`.

IDs must match abundance feature names. Confidence is `high`, `medium`, or `low`,
assigned by the analyst under an explicit evidence policy. These labels are not
calibrated probabilities. Combine multiple evidence sources for a pair in its
`evidence` field. A predicted host genus must not be presented as confirmed
species/strain specificity; link only at a resolution justified by the evidence.

## Configuration

```yaml
title: "Medication-associated microbiome trajectories"
synthetic: false
metadata: metadata.tsv
bacteria: bacteria.tsv
phages: phages.tsv          # optional; omit for bacteria-only studies
host_links: host_links.tsv  # optional; requires phages
output: results
abundance_scale: counts
phase_levels: [baseline, treatment, followup]
group_reference: comparator # omit for a single-cohort analysis
covariates: [age_baseline, bmi_baseline, metformin_baseline]
sensitivity_covariates: [weight_change_pct, fiber_change_g_day]
min_prevalence: 0.10
pseudocount: 0.000001
correlation: independent    # or car1; random intercept is always included
time_unit_days: 7           # CAR(1) correlation is per this many days
min_baseline_samples: 2
recovery_threshold: 0.15    # descriptive Bray-Curtis threshold, not validated
recovery_consecutive: 2
detection_threshold: 0     # relative-abundance units
host_confidence: [high, medium]
bootstrap_replicates: 500
seed: 20260928
report_engine: auto         # auto, quarto, or builtin
```

Covariates are a study-specific choice, not a universal prescription. Missing
required values, constant covariates, and exact collinearity are rejected. No
imputation or automated confounder selection is performed. The primary model uses
`covariates`; the separate sensitivity model also includes
`sensitivity_covariates`. Categorical reference levels follow R factor ordering
except for the explicitly configured phase and group reference levels.

Quarto is detected on PATH or in common RStudio locations. Set
`CHRONOMICROBIOME_QUARTO` to a full executable path if needed. In `auto` mode, an
unavailable or failing Quarto renderer produces an explicit fallback message and
the built-in standalone HTML report. The actual renderer is recorded in
`provenance.json`. `quarto` mode treats rendering failure as a failed run, even
though the built-in report is retained for inspection.

## Interpreting the model

For the two-group demo:

```text
CLR abundance ~ group * phase + age_baseline + bmi_baseline + metformin_baseline
Random intercept: participant_id
```

`groupsemaglutide:phasetreatment` is the semaglutide-minus-comparator difference
in baseline-to-treatment change, on the natural-log CLR scale. It is not a
clinical drug effect, an absolute abundance difference, or an ordinary fold change.
Without a comparator, phase coefficients cannot distinguish exposure-associated
changes from secular time changes.

See [Methods](docs/METHODS.md) for zero handling, multiple testing, model failures,
bootstrap limitations, and recovery definitions. See [Semaglutide study design](docs/SEMAGLUTIDE.md)
for metadata and estimand considerations.

## Tests

```r
install.packages("testthat")
```

```bash
R CMD build .
R CMD check --no-manual chronomicrobiome_0.1.0.tar.gz
```

The tests check sample alignment, invalid inputs, design identifiability, CLR
reference stability, injected signal direction, multiple-testing scope, optional
serial correlation, recovery logic, HTML escaping, provenance, and end-to-end
execution. Synthetic checks demonstrate software behavior; they do not establish
real-world false-positive control, clinical validity, or external reproducibility.

## Not implemented yet

- Raw FASTQ processing, assembly, or host prediction.
- Recurring exposures, crossover designs, multi-arm comparisons, or continuous dose-response models.
- Formal lagged interaction inference, mediation, causal inference, or clinical prediction.
- Absolute-abundance inference, time-varying confounding correction, or informative-dropout correction.
- A validated biological recovery threshold or externally benchmarked differential-abundance method.

Before use for a manuscript, define the estimand, review study design and
diagnostics, assess sensitivity to zero handling and model specification, and
benchmark against appropriate established microbiome methods and independent data.
