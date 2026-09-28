# Application: semaglutide and the human gut microbiome

Semaglutide (including the Ozempic brand) is a motivating medication exposure,
not a special statistical mode. The same framework can assess other medications,
antibiotic courses, and dietary interventions with suitable longitudinal data.

## A useful first research question

How do bacterial and phage relative-abundance patterns change after medication
initiation compared with an appropriate comparator, and how persistent are these
changes during the observed follow-up?

The software reports associations. A randomized comparison can strengthen an
effect interpretation when analyzed according to its design, but this software
alone does not establish causality. In observational studies, indication for
treatment, baseline disease, concurrent medication, and differential follow-up
can create confounding and selection bias.

## What to record

- Participant/sample identifiers, collection dates, index date, and body site.
- Generic medication and formulation, initiation, actual exposure at collection,
  recorded dose, adherence, interruptions, and discontinuation dates when known.
- Indication, baseline metabolic status, age, and baseline anthropometrics.
- Concurrent drugs, especially metformin and antibiotic exposure, with dates.
- Diet measurements, gastrointestinal symptoms, and weight measurements over time.
- Collection protocol, sequencing assay/batch, and sequencing depth.

No universal adjustment set is assumed. Define the question and causal structure
before selecting covariates. For example, post-initiation diet, weight, and
gastrointestinal changes may be mediators rather than baseline confounders.
Automatically adjusting for them can remove part of the association of interest
or introduce bias. ChronoµBiome keeps additional conditioning variables in a
separately labeled sensitivity model; it does not perform mediation analysis or
identify direct effects.

This version reads these exposures through analyst-defined phase labels and
covariates. It does not infer exposure windows from dates, estimate dose effects,
or implement marginal structural models for time-varying confounding. Dedicated
methods would be needed for recurrent exposure, adherence-adjusted causal
questions, and informative dropout.

## Synthetic demonstration

The included data contain 32 simulated participants and 320 samples with invented
bacterial/phage identities. Group labels are comparator and semaglutide. Artificial
abundance signals make it possible to verify that the workflow detects known
directions. The schedule includes an illustrative treatment window followed by
follow-up; it is not a recommended treatment or discontinuation schedule.

The primary model adjusts for simulated baseline age, BMI, and metformin status.
The sensitivity model additionally conditions on simulated weight and dietary
fiber changes. Neither the synthetic input nor resulting plots represent real
drug effects. Real microbiome changes should not be assumed to match this example.

## Recent research context

A [2026 preprint on semaglutide, dietary fiber, and post-cessation follow-up](https://www.medrxiv.org/content/10.64898/2026.05.02.26352300v1)
illustrates the relevance of the question. It is a preprint, not peer-reviewed
confirmation. Its results were not used to set the simulated biological effects.
The data are not bundled or analyzed here. Data access, participant metadata,
licensing, and suitability must be checked before any real-cohort benchmark.

A [registered 2026 semaglutide/metformin study](https://ensaiosclinicos.gov.br/rg/RBR-6pjdr2r)
also highlights the importance of distinguishing medication combinations. A trial
registration describes a planned study; it does not establish a microbiome effect.
