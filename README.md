# Hormesis and control-based no-effect estimates

This public research compendium contains the manuscript draft, exploratory design simulation and two case studies developed alongside the `bayesnec` Example 10 vignette. The manuscript proposes a post-peak no-significant-effect concentration (NSEC) for a hormetic concentration–response curve. The vignette in [`bayesnec`](https://github.com/open-AIMS/bayesnec) focuses on using the package.

## Contents

- `manuscript.qmd`: manuscript text, tables and figure code.
- `analysis/hormesis_simulation.R`: Gaussian plug-in design simulation. It does not test Bayesian NSEC coverage.
- `analysis/hormesis_negri.R`: local import of the Negri et al. marine workbook.
- `analysis/hormesis_case_fits.R`: candidate fits, sampler screening, targeted refits and exact leave-one-out scores for lettuce.
- `analysis/hormesis_fit_plots.R`: candidate and model-averaged curves, model weights, ECNSEC values and concentration-level checks from saved fits.
- `analysis/hormesis_toxval.R`: the NSEC in each direction (decreasing and increasing) from the saved fits, calculated with `toxval::nsec()`.
- `data/`: lettuce observations, simulation summaries, sampler diagnostics, posterior predictive checks at each tested concentration, fitted means at the tested concentrations and case-study estimates.
- `figures/`: candidate and model-averaged case-study figures generated from the saved fits.

## Reproduction

Run `analysis/hormesis_simulation.R` and then `analysis/hormesis_case_fits.R` from the repository root. The case-fit script imports the marine workbook, fits eleven candidate equations, screens the samplers, calculates exact leave-one-out scores where needed and calls `analysis/hormesis_fit_plots.R` to regenerate the figures and tables. `analysis/hormesis_toxval.R` then reads the saved fits and needs [`toxval`](https://github.com/open-AIMS/toxval) 2.0.0 or later, installed from GitHub. These Bayesian fits and refits can take considerable time. It requires `bayesnec`, `brms`, `rstan`, `ggplot2`, `loo`, `readxl` and their dependencies; Quarto and `knitr` render the manuscript. No package installation is performed by these scripts. The saved results used R 4.6.1, `brms` 2.23.0, `rstan` 2.32.7 and `bayesnec` fitting code at commit [`9b99f600`](https://github.com/open-AIMS/bayesnec/commit/9b99f600). The full provenance is in the manuscript's Reproducibility section.

The marine workbook is available from [AIMS](https://data.aims.gov.au/data-download/615e208d-f6c7-4b2f-b99e-d4fc50614e48/Negri_etal_2021_Experimental-data_Toxicity-thresholds-gas-condensate-oils.xlsx) under [CC BY-NC 3.0 Australia](https://researchdata.edu.au/condensate-toxicity-testing/1702929). Its individual observations and fitted objects are downloaded or generated into `ignore/` and are excluded from this repository. The marine figure contains attributed treatment means. The lettuce data were copied from the GPL-2 `drc` source, with provenance in the manuscript. Code in this compendium is GPL-2; the external workbook retains its own licence.
