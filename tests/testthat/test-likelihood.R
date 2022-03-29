test_that("negLogLik and gradient separate and parameters works", {
  pars <- exampleParams()
  # just a tweak of the parameters
  pars$gamma <- pars$gamma * 1.1
  pars$lambda_0 <- pars$lambda_0 * 0.9
  pars$theta <- pars$theta * 1.2
  R <- resSimAWX(n_thousands = 1, params = pars) # AWX dataframe
  l <- # log likelihood for these data and some random parameter values
    negLogLik(
      gamma = pars$gamma,
      lambda_0 = pars$lambda_0,
      theta = pars$theta,
      AWX = R
    )
  expect_length({
    l
  }, n = 1)
  expect_true({
    is.finite(l)
  }, TRUE)

  g_l0_t <-  GL0T(pars)
  ll <- negLogLik.vec(g_l0_t = g_l0_t,AWX = R)

  expect_identical(l,ll)

  gl <- gradLogLik(gamma = pars$gamma,
               lambda_0 = pars$lambda_0,
               theta = pars$theta,
               AWX = R)

  glv <- gradLogLik.vec(g_l0_t = GL0T(pars),AWX = R)

  expect_equal(gl,glv)
})

test_that("MLE are all positive",{
  pars <- exampleParams()
  R <- resSimAWX(n_thousands = 1, params = pars) # AWX dataframe

  expect_true(all( mleAll(R,params = pars) > 0 ))


})
