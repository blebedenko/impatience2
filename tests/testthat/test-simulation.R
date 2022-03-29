test_that("listParam returns all parameters", {
  expect_length(object = {
    pars <- listParams(gamma = 1, lambda_0 = 2,theta = 3,eta = 4,mu = 5,s = 6)
    pars
  },
  n = 6)
  expect_equal({
    pars <- listParams(gamma = 1, lambda_0 = 2,theta = 3,eta = 4,mu = 5,s = 6)
    sort(names(pars))
  },
  expected = c("eta", "gamma", "lambda_0", "mu", "s", "theta")) # sorted A-Z to prevent confusion

  expect_error(object = {listParams(1,1,1,1)}) # missing arguments = error

})



test_that("nextArrivalCosine gives correct time", {
  expect_true({nextArrivalCosine(current_time = 1.1,params = exampleParams()) > 1.1})

})


test_that("resSimCosine outputs a list", {
  expect_type({ resSimCosine(n = 5, params = exampleParams())},type = "list")
})


test_that("resSimAWX returns dataframe", {
  R <- resSimAWX(n_thousands = 1,params = exampleParams())
  expect_true({is.data.frame(R)})
  expect_length({names(R)}, n = 3)

})
