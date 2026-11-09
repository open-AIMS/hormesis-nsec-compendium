# Hormesis and control-based no-effect estimates

This public research compendium contains the manuscript draft, exploratory design simulation and two case studies developed alongside the `bayesnec` Example 10 vignette. The manuscript proposes a post-peak no-significant-effect concentration (NSEC) for a hormetic concentration–response curve. The vignette in [`bayesnec`](https://github.com/open-AIMS/bayesnec) focuses on using the package.

## Contents

- `manuscript.qmd`: manuscript text, tables and figure code.
- `analysis/hormesis_simulation.R`: Gaussian plug-in design simulation. It does not test Bayesian NSEC coverage.
- `analysis/hormesis_negri.R`: local import of the Negri et al. marine workbook.
- `analysis/hormesis_case_fits.R`: `bayesnec` fits and saved fitted objects.
- `analysis/hormesis_fit_plots.R`: case-study figures and concentration-level checks from saved fits.
- `data/`: lettuce observations and simulation summaries.
- `figures/`: case-study figures generated from the saved fits.

## Reproduction

Run scripts from the repository root in this order: `analysis/hormesis_simulation.R`, `analysis/hormesis_case_fits.R`, then `analysis/hormesis_fit_plots.R`. The case-fit script calls the local marine import and can take considerable time to fit eleven Bayesian models. It requires `bayesnec`, `brms`, `ggplot2`, `loo`, `readxl` and their dependencies; Quarto and `knitr` render the manuscript. No package installation is performed by these scripts. The version of R used for the saved fits was 4.6.1, with `brms` 2.23.0 and `rstan` 2.32.7. The full fit provenance is in the manuscript's Reproducibility section.

The marine workbook is available from [AIMS](https://data.aims.gov.au/data-download/615e208d-f6c7-4b2f-b99e-d4fc50614e48/Negri_etal_2021_Experimental-data_Toxicity-thresholds-gas-condensate-oils.xlsx) under [CC BY-NC 3.0 Australia](https://researchdata.edu.au/condensate-toxicity-testing/1702929). Its individual observations and fitted objects are downloaded or generated into `ignore/` and are excluded from this repository. The marine figure contains attributed treatment means. The lettuce data were copied from the GPL-2 `drc` source, with provenance in the manuscript. Code in this compendium is GPL-2; the external workbook retains its own licence.
