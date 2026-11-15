testthat::test_that("KM survival SE agrees with empirical binomial uncertainty", {
  d <- data.frame(time = 1:10, event = 1)
  ans <- summaryKM(d, "time", "event", int_name = "A",
                   types = "survival", se = TRUE)
  observed <- ans$value[ans$variable == "se" & ans$time == 8]
  testthat::expect_equal(observed, sqrt(0.2 * 0.8 / 10))
})

testthat::test_that("weighted survival SE respects both survfit uncertainty scales", {
  d <- data.frame(time = c(1, 2, 2, 4, 5, 6),
                  event = c(1, 0, 1, 1, 0, 1),
                  weight = c(0.5, 1.2, 0.8, 1.1, 0.7, 1.3))
  for (robust in c(FALSE, TRUE)) {
    sf <- survival::survfit(survival::Surv(time, event) ~ 1, data = d,
                           weights = weight, robust = robust)
    testthat::expect_identical(sf$logse, !robust)
    expected <- summary(sf, censored = TRUE)
    ans <- summaryKM(d, "time", "event", weight_var = "weight",
                     int_name = "A", types = "survival", se = TRUE,
                     robust = robust)
    observed <- ans[ans$variable == "se", ]
    testthat::expect_equal(observed$time, expected$time)
    testthat::expect_equal(observed$value, expected$std.err)
  }
})
