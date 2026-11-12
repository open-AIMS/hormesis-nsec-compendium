# Estimate the NSEC in each direction with toxval for the case-study fits, for
# each retained equation and for the model average. Run from the compendium
# root after analysis/hormesis_case_fits.R has saved the fits and model sets
# under ignore/hormesis_data/.

stopifnot(utils::packageVersion("toxval") >= "2.0.0")

fit_dir <- "ignore/hormesis_data"
sig_val <- 0.01

read_case <- function(case) {
  set <- readRDS(file.path(fit_dir, sprintf("%s_set_final.rds", case)))
  fits <- lapply(set$success_models, function(model) {
    readRDS(file.path(fit_dir, sprintf("fit_final_%s_%s.rds", case, model)))
  })
  c(list(`Model average` = set), stats::setNames(fits, set$success_models))
}

# toxval masks bayesnec::nsec(), so both are called with their prefixes. A
# realisation that never crosses its reference is Inf in value.draw.
summarise_fit <- function(fit) {
  two_way <- suppressWarnings(toxval::nsec(fit, sig_val = sig_val, draws = TRUE))
  post_peak <- suppressWarnings(bayesnec::nsec(fit, sig_val = sig_val))
  censored <- vapply(two_way$draws, function(d) {
    mean(is.infinite(d$value.draw[!is.na(d$value.draw)]))
  }, numeric(1))
  data.frame(
    method = c(paste("toxval", two_way$direction), "bayesnec post-peak"),
    estimate = c(two_way$estimate, post_peak[[1]]),
    conf.low = c(two_way$conf.low, post_peak[[2]]),
    conf.high = c(two_way$conf.high, post_peak[[3]]),
    censored = c(censored, NA),
    reference = c(two_way$reference, NA),
    control = c(two_way$control, NA)
  )
}

rows <- lapply(c("marine", "lettuce"), function(case) {
  members <- read_case(case)
  do.call(rbind, lapply(names(members), function(model) {
    cbind(case = case, model = model, summarise_fit(members[[model]]))
  }))
})
out <- do.call(rbind, rows)
print(out, digits = 3)
utils::write.csv(out, "data/hormesis_toxval_nsec.csv", row.names = FALSE)
