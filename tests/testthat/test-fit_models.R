test_that("fit model works as expected", {
  result <- fit_models()
  print("hi")
  expect_list(result,type="list")
})
