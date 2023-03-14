test_that("fit model works as expected", {
  result <- fit_models()
  expect_list(result,type="list")
})
