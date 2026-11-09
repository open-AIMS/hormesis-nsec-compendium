# Run from the compendium root after installing bayesnec and readxl.
library(bayesnec)
set.seed(570)
source("analysis/hormesis_negri.R")
lettuce <- read.csv("data/hormesis_lettuce_1993.csv")
negri <- read.csv("ignore/hormesis_data/hormesis_negri_2021.csv")
negri$conc <- (negri$source_conc_ug_l_tah - 1) / 1000
fit_negri <- bnec(sgr ~ crf(conc, model = "neclinhorme"),
                  data = negri, family = gaussian(),
                  chains = 4, cores = 4, iter = 2000, warmup = 1000,
                  seed = 570,
                  control = list(adapt_delta = 0.995, max_treedepth = 12))
curve_params(fit_negri)
nec(fit_negri) # the reported no-effect estimate
nsec(fit_negri, sig_val = 0.01)
ecx(fit_negri, ecx_val = 10)

fit_negri_smooth <- bnec(sgr ~ crf(conc, model = "ecxhormebc4"),
                         data = negri, family = gaussian(),
                         chains = 4, cores = 4, iter = 2000, warmup = 1000,
                         seed = 570, control = list(adapt_delta = 0.995,
                                                    max_treedepth = 12))
nsec(fit_negri_smooth, sig_val = 0.01)
ecx(fit_negri_smooth, ecx_val = 10)

fit_negri_monotone <- bnec(sgr ~ crf(conc, model = "neclin"),
                           data = negri, family = gaussian(),
                           chains = 4, cores = 4, iter = 2000, warmup = 1000,
                           seed = 570, control = list(adapt_delta = 0.995,
                                                      max_treedepth = 12))
nsec(fit_negri_monotone, sig_val = 0.01)
ecx(fit_negri_monotone, ecx_val = 10)
fit_marine_nechorme <- bnec(sgr ~ crf(conc, model = "nechorme"),
                            data = negri, family = gaussian(), chains = 4,
                            cores = 4, iter = 2000, warmup = 1000, seed = 570,
                            control = list(adapt_delta = 0.995,
                                           max_treedepth = 12))
fit_marine_nechorme4 <- bnec(sgr ~ crf(conc, model = "nechorme4"),
                             data = negri, family = gaussian(), chains = 4,
                             cores = 4, iter = 2000, warmup = 1000, seed = 570,
                             control = list(adapt_delta = 0.995,
                                            max_treedepth = 12))
fit_marine_bc5 <- bnec(sgr ~ crf(conc, model = "ecxhormebc5"),
                       data = negri, family = gaussian(), chains = 4,
                       cores = 4, iter = 2000, warmup = 1000, seed = 570,
                       init = 0, control = list(adapt_delta = 0.995,
                                                max_treedepth = 12))
marine_fits <- list(neclin = fit_negri_monotone,
                    neclinhorme = fit_negri,
                    nechorme = fit_marine_nechorme,
                    nechorme4 = fit_marine_nechorme4,
                    ecxhormebc4 = fit_negri_smooth,
                    ecxhormebc5 = fit_marine_bc5)
loo::loo_compare(list(
  neclin = loo::loo(fit_negri_monotone$fit),
  neclinhorme = loo::loo(fit_negri$fit),
  nechorme = loo::loo(fit_marine_nechorme$fit),
  nechorme4 = loo::loo(fit_marine_nechorme4$fit),
  ecxhormebc4 = loo::loo(fit_negri_smooth$fit),
  ecxhormebc5 = loo::loo(fit_marine_bc5$fit)
))

# Fitted and observed means at each tested concentration
treatments <- data.frame(conc = sort(unique(negri$conc)))
fitted_mean <- sapply(marine_fits, function(f) {
  apply(brms::posterior_epred(f$fit, newdata = treatments), 2, median)
})
cbind(treatments, observed = tapply(negri$sgr, negri$conc, mean), fitted_mean)

fit_lettuce <- bnec(weight ~ crf(log(conc + 0.1), model = "neclinhorme"),
                    data = lettuce, family = gaussian(),
                    chains = 4, cores = 4, iter = 2000, warmup = 1000,
                    seed = 570, control = list(adapt_delta = 0.99,
                                               max_treedepth = 12))
curve_params(fit_lettuce)
nsec(fit_lettuce, sig_val = 0.01)
ecx(fit_lettuce, ecx_val = 10)
fit_lettuce_nechorme <- bnec(weight ~ crf(log(conc + 0.1), model = "nechorme"),
                             data = lettuce, family = gaussian(), chains = 4,
                             cores = 4, iter = 2000, warmup = 1000, seed = 570,
                             control = list(adapt_delta = 0.99,
                                            max_treedepth = 12))
fit_lettuce_nechorme4 <- bnec(weight ~ crf(log(conc + 0.1), model = "nechorme4"),
                              data = lettuce, family = gaussian(), chains = 4,
                              cores = 4, iter = 2000, warmup = 1000, seed = 570,
                              control = list(adapt_delta = 0.99,
                                             max_treedepth = 12))
fit_lettuce_ecxhormebc4 <- bnec(weight ~ crf(log(conc + 0.1),
                                             model = "ecxhormebc4"),
                                data = lettuce, family = gaussian(), chains = 4,
                                cores = 4, iter = 2000, warmup = 1000,
                                seed = 570, control = list(adapt_delta = 0.99,
                                                           max_treedepth = 12))
fit_lettuce_bc5 <- bnec(weight ~ crf(log(conc + 0.1),
                                     model = "ecxhormebc5"),
                        data = lettuce, family = gaussian(), chains = 4,
                        cores = 4, iter = 2000, warmup = 1000, seed = 570,
                        init = 0, control = list(adapt_delta = 0.99,
                                                 max_treedepth = 12))
saveRDS(fit_negri, "ignore/hormesis_data/fit_neclinhorme.rds")
saveRDS(fit_negri_smooth, "ignore/hormesis_data/fit_ecxhormebc4.rds")
saveRDS(fit_negri_monotone, "ignore/hormesis_data/fit_neclin.rds")
saveRDS(fit_marine_nechorme, "ignore/hormesis_data/fit_marine_nechorme.rds")
saveRDS(fit_marine_nechorme4, "ignore/hormesis_data/fit_marine_nechorme4.rds")
saveRDS(fit_marine_bc5, "ignore/hormesis_data/fit_marine_bc5.rds")
saveRDS(fit_lettuce, "ignore/hormesis_data/fit_lettuce_log_neclinhorme.rds")
saveRDS(fit_lettuce_nechorme,
        "ignore/hormesis_data/fit_lettuce_log_nechorme.rds")
saveRDS(fit_lettuce_nechorme4,
        "ignore/hormesis_data/fit_lettuce_log_nechorme4.rds")
saveRDS(fit_lettuce_ecxhormebc4,
        "ignore/hormesis_data/fit_lettuce_log_ecxhormebc4.rds")
saveRDS(fit_lettuce_bc5,
        "ignore/hormesis_data/fit_lettuce_log_ecxhormebc5.rds")
source("analysis/hormesis_fit_plots.R")
