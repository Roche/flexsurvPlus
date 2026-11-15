test_that("model summaries report optimizer convergence without changing estimates", {
  data <- data.frame(time = c(1, 2, 3, 4, 6, 8, 10, 12),
                     event = c(1, 1, 0, 1, 1, 0, 1, 0))
  limited <- flexsurv::flexsurvreg(survival::Surv(time, event) ~ 1,
    data = data, dist = "weibull", inits = c(1, 5),
    control = list(maxit = 1))
  healthy <- flexsurv::flexsurvreg(survival::Surv(time, event) ~ 1,
    data = data, dist = "weibull")
  expect_identical(limited$opt$convergence, 1L)
  expect_identical(healthy$opt$convergence, 0L)
  models <- list(limited = limited, healthy = healthy)
  original <- models

  result <- get_model_summary(models)
  expect_identical(result$Status,
                   c("Not converged (optim code 1)", "Converged"))
  expect_identical(models, original)
  expect_equal(result$AIC, unname(vapply(models, stats::AIC, numeric(1))))
  expect_equal(result$BIC, unname(vapply(models, stats::BIC, numeric(1))))

  healthy$opt$convergence <- NULL
  expect_identical(get_model_summary(list(healthy = healthy))$Status,
                   "Unknown (optimizer code unavailable)")
})
