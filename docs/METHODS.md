# Methods and statistical scope

## Input contract

The program requires identical sample sets for metadata and supplied modalities,
unique sample identifiers and participant/time observations, finite nonnegative
abundances, and complete required covariates. It reorders abundance rows by IDs.
All-zero sample libraries are rejected. Counts must be integer; proportions must
sum to 1 within 0.01 and are renormalized exactly. Time is numeric in days.

At least six participants and two samples per participant are required. Two-group
comparisons require three participants per arm and observations in every
group/phase cell. These are engineering guardrails, **not power calculations or
proof of adequate inference**. Phases must progress in the configured order.
The first phase is baseline. Optional recovery uses a phase named followup.

## Transformation and filtering

For each modality separately, divide feature values by the sample total to obtain
relative abundance p. Remove globally all-zero features, add the configured
relative-unit pseudocount a, and transform:

`CLR_ij = log(p_ij + a) - mean_k(log(p_ik + a))`.

The reference mean uses all globally nonzero supplied features. Prevalence
filtering selects modeled features **after** constructing this reference, so
changing the prevalence cutoff alone does not change the CLR denominator.
Prevalence is the fraction of samples exceeding the detection threshold; sample
visits are not reweighted by participant. Unequal follow-up can influence this
filter. Zero replacement is a modeling assumption. Sensitivity to the pseudocount,
reference feature set, prevalence cutoff, and sequencing depth must be assessed.

Do not combine virus-enriched and bulk metagenomic abundance estimates as if they
had a common measurement scale. Both tables must represent comparable sampling
units and compatible assays. Separately normalized modalities cannot establish
absolute phage/host abundance ratios.

## Exposure-associated models

`nlme::lme` fits REML linear mixed models with a participant random intercept.
Fixed effects are phase, configured covariates, and—when a comparator is
specified—group and group-by-phase interaction. Time is used for plotting,
ordering, and optional residual correlation; the fixed exposure effect is
**phase-based**, not a continuously varying trajectory smoother.

`correlation: independent` assumes conditionally independent residuals after the
random intercept. `car1` uses `nlme::corCAR1` with actual time divided by
`time_unit_days`, allowing irregular visit spacing. The latter assumes a common
continuous-time AR(1) structure, not arbitrary autocorrelation.

Two-group primary contrasts are group-by-nonbaseline-phase interactions. For a
single cohort, they are phase coefficients. Both use the first configured phase
as baseline. Causal difference-in-differences interpretation would require
additional design assumptions, including suitable counterfactual trends; the
program does not test or guarantee them.

Pointwise 95% intervals use the coefficient-specific degrees of freedom returned
by nlme. BH correction is applied jointly across all planned primary
feature-by-phase tests in bacteria and phages, separately for primary and
sensitivity specifications. The planned family includes failed model fits.
Other coefficients are retained with raw p-values but no q-values.

No successful fit is fabricated when a model fails. Failed/nonvarying features
remain in diagnostics. Convergence warnings and near-boundary random-intercept
variance are reported. If every model fails, the run fails. When only some fits
fail, the report retains successful results and the diagnostics show failures;
completion is not a guarantee of model adequacy. A user must inspect residuals,
heteroscedasticity, influential participants, and serial dependence. Provided
diagnostic plots cover a small deterministic subset; a full residual table is saved.

This is an exploratory CLR association workflow, not a new validated
differential-abundance estimator. In particular, sparse data can violate the
Gaussian working model, and BH adjustment does not repair misspecified p-values.

## Persistence and recovery

Persistence is an observed detection fraction per participant/feature, with first
and last detected sampling dates. It does not imply persistence between visits.

The personal baseline reference is the mean relative composition of baseline
samples. Each baseline sample is excluded from its own reference. Bray–Curtis
distance between normalized compositions is half their L1 distance.

Recovery is descriptive and requires:

1. At least `min_baseline_samples` baseline samples.
2. At least one observation in an exposure phase (neither baseline nor followup).
3. An observed exposure-phase distance strictly above `recovery_threshold`.
4. `recovery_consecutive` successive observed follow-up samples at or below that
   threshold.

The returned date is the first sample of the qualifying sequence; confirmation
is the last sample in that sequence. If available, the previous follow-up sample
date is retained, but no exact event time or continuous-time recovery curve is
estimated. Long gaps between qualifying samples are not bounded in v0.1. Users
must judge whether their sampling density supports interpretation.

Participants without an observed departure, without follow-up, with inadequate
baseline, or without a qualifying sequence receive different statuses. No
observed return by the last sample does not imply permanent failure to recover.
The default threshold is an illustrative analyst setting, not a biological norm.

## Phage–host associations

Only supplied links are evaluated. For each pair, calculate a Spearman correlation
of the two primary models' normalized residuals within each participant. Require
four observations and variable residuals. Average these correlations with equal
participant weights; require six eligible participants for a summary. Resample
participant correlations to obtain percentile intervals.

Intervals condition on fitted residuals: models are not refitted during the
bootstrap. Intervals omit model-fitting and host-assignment uncertainty and are
not multiplicity-adjusted. There is no significance call or p-value. Participants
from both groups are pooled in this descriptive summary; treatment-specific
changes in association are not tested. Shared compositional artifacts may induce
correlations. These outputs do not demonstrate infection, killing, causation, or
temporal direction.

Confidence sensitivity displays the configured primary evidence subset and the
high-only subset. It does not average prediction scores or model host uncertainty
probabilistically. The analyst must define and document the evidence categories.

## Provenance and reporting

Inputs receive MD5 reproducibility fingerprints (not cryptographic security
guarantees). Resolved configuration, UTC times, software versions, RNG seed,
model failures, tables, and session information are saved. Seeded simulation and
bootstrap calls preserve the caller's R random-number state. Quarto renders when
available; the built-in fallback uses the same report content. No report assets
are fetched from the network.

## References

- [nlme mixed models](https://stat.ethz.ch/R-manual/R-devel/library/nlme/html/lme.html)
- [Continuous-time AR(1)](https://stat.ethz.ch/R-manual/R-devel/library/nlme/html/corCAR1.html)
- [R multiple-comparison adjustments](https://stat.ethz.ch/R-manual/R-devel/library/stats/html/p.adjust.html)
- [MaAsLin 3](https://huttenhower.sph.harvard.edu/maaslin) is an established related
  method to consider in external benchmarking; it is not integrated in v0.1.

Research references and naming checks were reviewed on 2026-09-28.
