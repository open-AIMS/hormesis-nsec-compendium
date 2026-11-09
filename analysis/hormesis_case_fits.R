# Run from the compendium root after installing bayesnec and readxl.
library(bayesnec)
options(brms.backend = "rstan", mc.cores = 4)
rstan::rstan_options(auto_write = TRUE)
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

fit_negri_smooth <- bnec(sgr ~ crf(conc, model = "ecxhormebc4"),
                         data = negri, family = gaussian(),
                         chains = 4, cores = 4, iter = 2000, warmup = 1000,
                         seed = 570, control = list(adapt_delta = 0.995,
                                                    max_treedepth = 12))

fit_negri_monotone <- bnec(sgr ~ crf(conc, model = "neclin"),
                           data = negri, family = gaussian(),
                           chains = 4, cores = 4, iter = 2000, warmup = 1000,
                           seed = 570, control = list(adapt_delta = 0.995,
                                                      max_treedepth = 12))
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
fit_lettuce <- bnec(weight ~ crf(log(conc + 0.1), model = "neclinhorme"),
                    data = lettuce, family = gaussian(),
                    chains = 4, cores = 4, iter = 2000, warmup = 1000,
                    seed = 570, control = list(adapt_delta = 0.99,
                                               max_treedepth = 12))
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
                        init = 0, control = list(adapt_delta = 0.999,
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

# Retain equations whose final sampler has Rhat <= 1.01, ESS >= 400 and no
# divergent transitions. Longer runs preserve 4000 post-warm-up draws per fit.
refit_longer <- function(fit, adapt_delta) {
  update(fit, iter = 6000, warmup = 4000, thin = 2, seed = 570,
         control = list(adapt_delta = adapt_delta, max_treedepth = 15),
         refresh = 0)
}
fit_negri <- refit_longer(fit_negri, 0.999)
fit_marine_bc5 <- refit_longer(fit_marine_bc5, 0.999)
fit_lettuce_nechorme4 <- refit_longer(fit_lettuce_nechorme4, 0.999)
fit_lettuce_ecxhormebc4 <- refit_longer(fit_lettuce_ecxhormebc4, 0.99)
fit_lettuce_ecxhormebc4 <- refit_longer(fit_lettuce_ecxhormebc4, 0.9999)
fit_lettuce_bc5 <- refit_longer(fit_lettuce_bc5, 0.9999)

marine_fits <- list(neclin = fit_negri_monotone,
                    neclinhorme = fit_negri,
                    nechorme = fit_marine_nechorme,
                    nechorme4 = fit_marine_nechorme4,
                    ecxhormebc4 = fit_negri_smooth,
                    ecxhormebc5 = fit_marine_bc5)
lettuce_fits <- list(neclinhorme = fit_lettuce,
                     nechorme = fit_lettuce_nechorme,
                     nechorme4 = fit_lettuce_nechorme4,
                     ecxhormebc4 = fit_lettuce_ecxhormebc4,
                     ecxhormebc5 = fit_lettuce_bc5)
retain_sampled <- function(fits, case) {
  diagnostics <- do.call(rbind, lapply(names(fits), function(model) {
    data.frame(case = case,
               check_sampling(fits[[model]], divergence_cutoff = 0))
  }))
  kept <- fits[!diagnostics$failed]
  if (length(kept) < 2) stop("Fewer than two equations passed the sampler screen")
  set.seed(570)
  set <- screen_models(do.call(base::c, unname(kept)),
                       divergence_cutoff = 0)
  list(set = set, diagnostics = diagnostics)
}
marine_result <- retain_sampled(marine_fits, "marine")
lettuce_result <- retain_sampled(lettuce_fits, "lettuce")
marine_set <- marine_result$set
set.seed(570)
lettuce_set <- amend(lettuce_result$set,
                     loo_controls = list(fitting = list(reloo = TRUE)))
stopifnot(all(!check_sampling(marine_set, divergence_cutoff = 0)$failed),
          all(!check_sampling(lettuce_set, divergence_cutoff = 0)$failed))
write.csv(rbind(marine_result$diagnostics, lettuce_result$diagnostics),
          "data/hormesis_sampler_diagnostics.csv", row.names = FALSE)
saveRDS(marine_set, "ignore/hormesis_data/marine_set_final.rds")
saveRDS(lettuce_set, "ignore/hormesis_data/lettuce_set_final.rds")
for (model in names(marine_fits)[!marine_result$diagnostics$failed]) {
  saveRDS(marine_fits[[model]], file.path("ignore/hormesis_data",
                                        paste0("fit_final_marine_", model,
                                               ".rds")))
}
for (model in names(lettuce_fits)[!lettuce_result$diagnostics$failed]) {
  saveRDS(lettuce_fits[[model]], file.path("ignore/hormesis_data",
                                         paste0("fit_final_lettuce_", model,
                                                ".rds")))
}
source("analysis/hormesis_fit_plots.R")
