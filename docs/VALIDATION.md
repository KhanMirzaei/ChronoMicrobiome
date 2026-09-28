# Validation — ChronoµBiome 0.1.0

Validated locally on 2026-09-28 with R 4.5.2 on macOS arm64.

- R CMD build: passed.
- R CMD check --no-manual: Status: OK (0 errors, 0 warnings, 0 notes in final status).
- Automated tests: 50 passing assertions, 0 failures, 0 warnings, 0 skips.
- Source CLI doctor, demo, validate, and run: exercised.
- Synthetic demo: 32 participants, 320 samples, 12 bacterial features, 8 phage features.
- All 40 primary/sensitivity mixed models fitted. A near-boundary random intercept
  for sparse Bacterium_12 is explicitly flagged in both specifications.
- Quarto rendering: passed using bundled Quarto after allowing its standard local cache.
- Built-in standalone HTML fallback: passed in the end-to-end automated test.
- Report structure: six embedded images, all required sections, no broken local download links.
- Trajectory and contrast figures were visually inspected. Full browser-layout QA
  was unavailable because the browser tool blocked local file URLs; no workaround
  was used. Browser appearance across devices has not been certified.

The check environment could not refresh CRAN/Bioconductor repository indexes;
checks used installed dependencies and still completed with Status: OK.

This release has NOT been validated on a real cohort, externally benchmarked for
statistical calibration, or reviewed for clinical use. The synthetic semaglutide
example is not evidence about drug effects. The display name is ChronoµBiome and the package is chronomicrobiome.
An unrelated project uses Chronobiome (chronobiome.org); this project is independent.
