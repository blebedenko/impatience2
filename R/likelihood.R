#' (negative) mean log-likelihood for the sinusoidal arrivals model
#' The mean is returned instead of the sum - should help gradient based optimizers
#' @param gamma periodical component of rate function
#' @param lambda_0 constant component of rate function
#' @param theta parameter of exponential patience
#' @param AWX compact simulation results for memory economy (A,W,X).
#' @return The negative log-likelihood at the point provided.
#' @export
negLogLik <- function(gamma,lambda_0,theta, AWX) {
  # gamma <- params$gamma
  # lambda_0 <- params$lambda_0
  # theta <- params$theta

  A <- AWX$A
  W <- AWX$W
  X <- AWX$X
  A_i = A[-1]
  A_tilde <- c(0, cumsum(A))
  A_tilde <- A_tilde[-length(A_tilde)]
  A_tilde_i = cumsum(A_i)
  W_i = W[-1]
  w_i = W[-length(W)]
  x_i = X[-length(X)]

  # elements of the log-likelihood
  l_i <-
    log(gamma / 2 + lambda_0 + (gamma * cos(A_tilde_i * pi * 2)) / 2) +
    log(exp(-W_i * theta)) +
    (gamma * exp(-theta * (w_i + x_i)) * (
      pi * sin(A_tilde_i * pi * 2) * 2 + theta * cos(A_tilde_i * pi * 2) - pi *
        sin(pi * (A_i + A_tilde_i) * 2) * exp(A_i * theta) * 2 - theta * exp(A_i *
                                                                               theta) * cos(pi * (A_i + A_tilde_i) * 2)
    )) / (pi ^ 2 * 8 + theta ^ 2 * 2) - (lambda_0 * exp(-theta * (w_i + x_i)) *
                                           (exp(A_i * theta) - 1)) / theta - (gamma * exp(-theta * (w_i + x_i)) * (exp(A_i *
                                                                                                                         theta) - 1)) / (theta * 2)

  # return the negative mean
  negMean <- -mean(l_i)
  return(negMean)

}


#' Wrapper for negLogLik to use with optim
#'
#' @param g_l0_t vector c(gamma, lambda_0, theta)
#' @param AWX AWX dataframe
#' @param which_known if supplied, supposedly returns a partial likelihood where this is fixed.
#' @return negative mean likelihood for the point g_l0_t
negLogLik.vec <- function(g_l0_t, AWX, which_known = NULL){

  gamma <- g_l0_t[1]
  lambda_0 <- g_l0_t[2]
  theta <- g_l0_t[3]
  if (is.null(which_known))
    return(negLogLik(gamma = gamma, lambda_0 = lambda_0 , theta = theta, AWX = AWX))
  # one known parameter
  # if (length(which_known == 1)){
  #
  # }

}


#' Gradient of the full likelihood, or a partial version
#'
#' @param gamma gamma
#' @param lambda_0 lambda_0
#' @param theta theta
#' @param AWX AWX
#'
#' @return the gradient at (gamma,lambda_0,theta) for data AWX
#' @export
gradLogLik <- function(gamma,lambda_0,theta, AWX) {

  A <- AWX$A
  W <- AWX$W
  X <- AWX$X
  A_i = A[-1]
  A_tilde <- c(0, cumsum(A))
  A_tilde <- A_tilde[-length(A_tilde)]
  A_tilde_i = cumsum(A_i)
  W_i = W[-1]
  w_i = W[-length(W)]
  x_i = X[-length(X)]



  # derivative by gamma:
  dl_gamma <-
    (cos(A_tilde_i * pi * 2) / 2 + 1 / 2) /
    (gamma / 2 + lambda_0 + (gamma * cos(A_tilde_i * pi * 2)) / 2) - (exp(-theta *
                                                                            (w_i + x_i)) * (exp(A_i * theta) - 1)) / (theta * 2) + (exp(-theta * (w_i +
                                                                                                                                                    x_i)) * (
                                                                                                                                                      pi * sin(A_tilde_i * pi * 2) * 2 + theta * cos(A_tilde_i * pi * 2) - pi *
                                                                                                                                                        sin(pi * (A_i + A_tilde_i) * 2) * exp(A_i * theta) * 2 - theta * exp(A_i *
                                                                                                                                                                                                                               theta) * cos(pi * (A_i + A_tilde_i) * 2)
                                                                                                                                                    )) / (pi ^ 2 * 8 + theta ^ 2 * 2)

  #derivative by lambda_0:
  dl_lambda_0 <-
    1 / (gamma / 2 + lambda_0 + (gamma * cos(A_tilde_i * pi * 2)) / 2) -
    (exp(-theta * (w_i + x_i)) * (exp(A_i * theta) - 1)) / theta

  # derivative by theta:
  dl_theta <-
    -W_i + lambda_0 * 1 / theta ^ 2 * exp(-theta * (w_i + x_i)) * (exp(A_i *
                                                                         theta) - 1) - (gamma * exp(-theta * (w_i + x_i)) * (
                                                                           -cos(A_tilde_i * pi * 2) + exp(A_i * theta) * cos(pi * (A_i + A_tilde_i) *
                                                                                                                               2) + A_i * pi * sin(pi * (A_i + A_tilde_i) * 2) * exp(A_i * theta) * 2 +
                                                                             A_i * theta * exp(A_i * theta) * cos(pi * (A_i + A_tilde_i) * 2)
                                                                         )) / (pi ^ 2 * 8 + theta ^ 2 * 2) + (gamma * 1 / theta ^ 2 * exp(-theta *
                                                                                                                                            (w_i + x_i)) * (exp(A_i * theta) - 1)) / 2 + (gamma * exp(-theta * (w_i +
                                                                                                                                                                                                                  x_i)) * (w_i + x_i) * (exp(A_i * theta) - 1)) / (theta * 2) - (
                                                                                                                                                                                                                    gamma * exp(-theta * (w_i + x_i)) * (w_i + x_i) * (
                                                                                                                                                                                                                      pi * sin(A_tilde_i * pi * 2) * 2 + theta * cos(A_tilde_i * pi * 2) - pi *
                                                                                                                                                                                                                        sin(pi * (A_i + A_tilde_i) * 2) * exp(A_i * theta) * 2 - theta * exp(A_i *
                                                                                                                                                                                                                                                                                               theta) * cos(pi * (A_i + A_tilde_i) * 2)
                                                                                                                                                                                                                    )
                                                                                                                                                                                                                  ) / (pi ^ 2 * 8 + theta ^ 2 * 2) + (lambda_0 * exp(-theta * (w_i + x_i)) *
                                                                                                                                                                                                                                                        (w_i + x_i) * (exp(A_i * theta) - 1)) / theta - gamma * theta * exp(-theta *
                                                                                                                                                                                                                                                                                                                              (w_i + x_i)) * 1 / (pi ^ 2 * 4 + theta ^ 2) ^ 2 * (
                                                                                                                                                                                                                                                                                                                                pi * sin(A_tilde_i * pi * 2) * 2 + theta * cos(A_tilde_i * pi * 2) - pi *
                                                                                                                                                                                                                                                                                                                                  sin(pi * (A_i + A_tilde_i) * 2) * exp(A_i * theta) * 2 - theta * exp(A_i *
                                                                                                                                                                                                                                                                                                                                                                                                         theta) * cos(pi * (A_i + A_tilde_i) * 2)
                                                                                                                                                                                                                                                                                                                              ) - (A_i * gamma * exp(A_i * theta) * exp(-theta * (w_i + x_i))) / (theta *
                                                                                                                                                                                                                                                                                                                                                                                                    2) - (A_i * lambda_0 * exp(A_i * theta) * exp(-theta * (w_i + x_i))) / theta
  # return the negative of the gradient elements' mean
  negativeGradientMean <-
    -c(mean(dl_gamma), mean(dl_lambda_0), mean(dl_theta))

  names(negativeGradientMean) <- c("gamma","lambda_0","theta")

  return(negativeGradientMean)




}

#' Gradient of  likelihood, vector form
#'
#' @param g_l0_t vector c(gamma, lambda_0, theta)
#' @param AWX AWX
#'
#' @return the gradient at (gamma,lambda_0,theta) for data AWX
#' @export
gradLogLik.vec <- function(g_l0_t, AWX) {
  gamma <- g_l0_t[1]
  lambda_0 <- g_l0_t[2]
  theta <- g_l0_t[3]
  A <- AWX$A
  W <- AWX$W
  X <- AWX$X
  A_i = A[-1]
  A_tilde <- c(0, cumsum(A))
  A_tilde <- A_tilde[-length(A_tilde)]
  A_tilde_i = cumsum(A_i)
  W_i = W[-1]
  w_i = W[-length(W)]
  x_i = X[-length(X)]



  # derivative by gamma:
  dl_gamma <-
    (cos(A_tilde_i * pi * 2) / 2 + 1 / 2) /
    (gamma / 2 + lambda_0 + (gamma * cos(A_tilde_i * pi * 2)) / 2) - (exp(-theta *
                                                                            (w_i + x_i)) * (exp(A_i * theta) - 1)) / (theta * 2) + (exp(-theta * (w_i +
                                                                                                                                                    x_i)) * (
                                                                                                                                                      pi * sin(A_tilde_i * pi * 2) * 2 + theta * cos(A_tilde_i * pi * 2) - pi *
                                                                                                                                                        sin(pi * (A_i + A_tilde_i) * 2) * exp(A_i * theta) * 2 - theta * exp(A_i *
                                                                                                                                                                                                                               theta) * cos(pi * (A_i + A_tilde_i) * 2)
                                                                                                                                                    )) / (pi ^ 2 * 8 + theta ^ 2 * 2)

  #derivative by lambda_0:
  dl_lambda_0 <-
    1 / (gamma / 2 + lambda_0 + (gamma * cos(A_tilde_i * pi * 2)) / 2) -
    (exp(-theta * (w_i + x_i)) * (exp(A_i * theta) - 1)) / theta

  # derivative by theta:
  dl_theta <-
    -W_i + lambda_0 * 1 / theta ^ 2 * exp(-theta * (w_i + x_i)) * (exp(A_i *
                                                                         theta) - 1) - (gamma * exp(-theta * (w_i + x_i)) * (
                                                                           -cos(A_tilde_i * pi * 2) + exp(A_i * theta) * cos(pi * (A_i + A_tilde_i) *
                                                                                                                               2) + A_i * pi * sin(pi * (A_i + A_tilde_i) * 2) * exp(A_i * theta) * 2 +
                                                                             A_i * theta * exp(A_i * theta) * cos(pi * (A_i + A_tilde_i) * 2)
                                                                         )) / (pi ^ 2 * 8 + theta ^ 2 * 2) + (gamma * 1 / theta ^ 2 * exp(-theta *
                                                                                                                                            (w_i + x_i)) * (exp(A_i * theta) - 1)) / 2 + (gamma * exp(-theta * (w_i +
                                                                                                                                                                                                                  x_i)) * (w_i + x_i) * (exp(A_i * theta) - 1)) / (theta * 2) - (
                                                                                                                                                                                                                    gamma * exp(-theta * (w_i + x_i)) * (w_i + x_i) * (
                                                                                                                                                                                                                      pi * sin(A_tilde_i * pi * 2) * 2 + theta * cos(A_tilde_i * pi * 2) - pi *
                                                                                                                                                                                                                        sin(pi * (A_i + A_tilde_i) * 2) * exp(A_i * theta) * 2 - theta * exp(A_i *
                                                                                                                                                                                                                                                                                               theta) * cos(pi * (A_i + A_tilde_i) * 2)
                                                                                                                                                                                                                    )
                                                                                                                                                                                                                  ) / (pi ^ 2 * 8 + theta ^ 2 * 2) + (lambda_0 * exp(-theta * (w_i + x_i)) *
                                                                                                                                                                                                                                                        (w_i + x_i) * (exp(A_i * theta) - 1)) / theta - gamma * theta * exp(-theta *
                                                                                                                                                                                                                                                                                                                              (w_i + x_i)) * 1 / (pi ^ 2 * 4 + theta ^ 2) ^ 2 * (
                                                                                                                                                                                                                                                                                                                                pi * sin(A_tilde_i * pi * 2) * 2 + theta * cos(A_tilde_i * pi * 2) - pi *
                                                                                                                                                                                                                                                                                                                                  sin(pi * (A_i + A_tilde_i) * 2) * exp(A_i * theta) * 2 - theta * exp(A_i *
                                                                                                                                                                                                                                                                                                                                                                                                         theta) * cos(pi * (A_i + A_tilde_i) * 2)
                                                                                                                                                                                                                                                                                                                              ) - (A_i * gamma * exp(A_i * theta) * exp(-theta * (w_i + x_i))) / (theta *
                                                                                                                                                                                                                                                                                                                                                                                                    2) - (A_i * lambda_0 * exp(A_i * theta) * exp(-theta * (w_i + x_i))) / theta
  # return the negative of the gradient elements' mean
  negativeGradientMean <-
    -c(mean(dl_gamma), mean(dl_lambda_0), mean(dl_theta))

  names(negativeGradientMean) <- c("gamma","lambda_0","theta")

  return(negativeGradientMean)




}
#' Cheap evaluation of the gradient
#'
#' @param params list of parameters
#' @param which_known vector of length 1-3, 1 = d_gamma 2 = d_lambda_0 3 = d_theta
#' @param AWX data
#'
#' @return a named vector with the corresponding gradient for the unknown parameters
grad.cheap <- function(params,which_known,AWX){

  gamma <- params$gamma
  lambda_0 <- params$lambda_0
  theta <- params$theta

  A <- AWX$A
  W <- AWX$W
  X <- AWX$X
  A_i = A[-1]
  A_tilde <- c(0, cumsum(A))
  A_tilde <- A_tilde[-length(A_tilde)]
  A_tilde_i = cumsum(A_i)
  W_i = W[-1]
  w_i = W[-length(W)]
  x_i = X[-length(X)]

  ANS <- list()
  # if gamma is unknown:
  if (!("gamma" %in% which_known)){

    # derivative by gamma:
    dl_gamma <-
      (cos(A_tilde_i * pi * 2) / 2 + 1 / 2) /
      (gamma / 2 + lambda_0 + (gamma * cos(A_tilde_i * pi * 2)) / 2) - (exp(-theta *
                                                                              (w_i + x_i)) * (exp(A_i * theta) - 1)) / (theta * 2) + (exp(-theta * (w_i +
                                                                                                                                                      x_i)) * (
                                                                                                                                                        pi * sin(A_tilde_i * pi * 2) * 2 + theta * cos(A_tilde_i * pi * 2) - pi *
                                                                                                                                                          sin(pi * (A_i + A_tilde_i) * 2) * exp(A_i * theta) * 2 - theta * exp(A_i *
                                                                                                                                                                                                                                 theta) * cos(pi * (A_i + A_tilde_i) * 2)
                                                                                                                                                      )) / (pi ^ 2 * 8 + theta ^ 2 * 2)

  ANS$dl_gamma <- dl_gamma
    }
    # if lambda_0 is unknown:
    if (!("lambda_0" %in% which_known)){
    #derivative by lambda_0:
    dl_lambda_0 <-
      1 / (gamma / 2 + lambda_0 + (gamma * cos(A_tilde_i * pi * 2)) / 2) -
      (exp(-theta * (w_i + x_i)) * (exp(A_i * theta) - 1)) / theta

    ANS$dl_lambda_0 <- dl_lambda_0

  }
  # if theta is unknown:
    if (!("theta" %in% which_known)){
    # derivative by theta:
    dl_theta <-
      -W_i + lambda_0 * 1 / theta ^ 2 * exp(-theta * (w_i + x_i)) * (exp(A_i *
                                                                           theta) - 1) - (gamma * exp(-theta * (w_i + x_i)) * (
                                                                             -cos(A_tilde_i * pi * 2) + exp(A_i * theta) * cos(pi * (A_i + A_tilde_i) *
                                                                                                                                 2) + A_i * pi * sin(pi * (A_i + A_tilde_i) * 2) * exp(A_i * theta) * 2 +
                                                                               A_i * theta * exp(A_i * theta) * cos(pi * (A_i + A_tilde_i) * 2)
                                                                           )) / (pi ^ 2 * 8 + theta ^ 2 * 2) + (gamma * 1 / theta ^ 2 * exp(-theta *
                                                                                                                                              (w_i + x_i)) * (exp(A_i * theta) - 1)) / 2 + (gamma * exp(-theta * (w_i +
                                                                                                                                                                                                                    x_i)) * (w_i + x_i) * (exp(A_i * theta) - 1)) / (theta * 2) - (
                                                                                                                                                                                                                      gamma * exp(-theta * (w_i + x_i)) * (w_i + x_i) * (
                                                                                                                                                                                                                        pi * sin(A_tilde_i * pi * 2) * 2 + theta * cos(A_tilde_i * pi * 2) - pi *
                                                                                                                                                                                                                          sin(pi * (A_i + A_tilde_i) * 2) * exp(A_i * theta) * 2 - theta * exp(A_i *
                                                                                                                                                                                                                                                                                                 theta) * cos(pi * (A_i + A_tilde_i) * 2)
                                                                                                                                                                                                                      )
                                                                                                                                                                                                                    ) / (pi ^ 2 * 8 + theta ^ 2 * 2) + (lambda_0 * exp(-theta * (w_i + x_i)) *
                                                                                                                                                                                                                                                          (w_i + x_i) * (exp(A_i * theta) - 1)) / theta - gamma * theta * exp(-theta *
                                                                                                                                                                                                                                                                                                                                (w_i + x_i)) * 1 / (pi ^ 2 * 4 + theta ^ 2) ^ 2 * (
                                                                                                                                                                                                                                                                                                                                  pi * sin(A_tilde_i * pi * 2) * 2 + theta * cos(A_tilde_i * pi * 2) - pi *
                                                                                                                                                                                                                                                                                                                                    sin(pi * (A_i + A_tilde_i) * 2) * exp(A_i * theta) * 2 - theta * exp(A_i *
                                                                                                                                                                                                                                                                                                                                                                                                           theta) * cos(pi * (A_i + A_tilde_i) * 2)
                                                                                                                                                                                                                                                                                                                                ) - (A_i * gamma * exp(A_i * theta) * exp(-theta * (w_i + x_i))) / (theta *
                                                                                                                                                                                                                                                                                                                                                                                                      2) - (A_i * lambda_0 * exp(A_i * theta) * exp(-theta * (w_i + x_i))) / theta
    ANS$dl_theta <- dl_theta

    }

  grad_vec <- Reduce(mean,ANS)
  return(grad_vec * (-1))

}


#' Function factory for partial likelihoods
#'
#' @param which_known Which of the parameter values is _known_ to the estimator? "theta", "gamma" or "lambda_0".
#' @param params vector c(gamma, lambda_0, theta)
#' @param AWX data
#'
#' @return a function that computes likelihood for the remaining two parameter except which_known.
#' The function's input is always a vector where the parameters are ordered alphabetically by their names, gamma before lambda_0 before theta.
#' @export
#'
#' @examples
#' params <- exampleParams()
#' AWX <- resSimAWX(4,params)
#' f <- oneKnownLik("theta",params = params,AWX = AWX)
#' f(c(1,2)) # likelihood for gamma = 1 lambda_0 = 2
oneKnownLik <- function(which_known,params, AWX){
  gamma <- params$gamma
  lambda_0 <- params$lambda_0
  theta <- params$theta

  if (which_known == "gamma") {
    negLik <- Vectorize(
      FUN = negLogLik,
      vectorize.args = c("lambda_0", "theta"),
      SIMPLIFY = TRUE
    )
    # likelihood for lambda_0, theta with known gamma
    negLik <- purrr::partial(negLik, AWX = AWX, gamma = gamma)
    oneLik <- function(lambda_0,theta) negLik(lambda_0,theta)

  }
  if (which_known == "lambda_0") {
    negLik <- Vectorize(
      FUN = negLogLik,
      vectorize.args = c("gamma", "theta"),
      SIMPLIFY = TRUE
    )
    # likelihood for gamma, theta with known lambda_0
    negLik <- purrr::partial(negLik, AWX = AWX, lambda_0 = lambda_0)
    oneLik <- function(gamma,theta) negLik(gamma,theta)

  }
  if (which_known == "theta") {
    negLik <- Vectorize(
      FUN = negLogLik,
      vectorize.args = c("gamma", "lambda_0"),
      SIMPLIFY = TRUE
    )
    # likelihood for gamma, lambda_0 with known theta
    negLik <- purrr::partial(negLik, AWX = AWX, theta = theta)
    oneLik <- function(gamma,lambda_0) negLik(gamma,lambda_0)
  }

  curr_params <- setdiff(c("gamma","lambda_0","theta"), which_known)

  negLogL <- function(two_par_vec) oneLik(two_par_vec[1],two_par_vec[2]) # wrap it for optim
  return(negLogL)
}

#' Title f
#'
#' @param which_known name of known parameter
#' @param params list of parameters
#' @param AWX dataset
#'
#' @return a gradient function for the two other parameter
#' @export
oneKnownGrad <- function(which_known,params, AWX){
  # gamma <- params$gamma
  # lambda_0 <- params$lambda_0
  # theta <- params$theta
  which_uknown <- setdiff(names(GL0T(params)), which_known)
  if (which_known == "gamma") {

    # gradient for lambda_0, theta with known gamma
    oneGrad <- purrr::partial(gradLogLik, AWX = AWX, gamma = params$gamma)

    }

  if (which_known == "lambda_0") {
    oneGrad <- purrr::partial(gradLogLik, AWX = AWX, lambda_0 = params$lambda_0)

  }
  if (which_known == "theta") {

    oneGrad <- purrr::partial(gradLogLik, AWX = AWX, theta = params$theta)

  }


  return(oneGrad)
}



#' Function factory for partial likelihoods
#'
#' @param which_not_known name of the UNKNOWN parameter to be estimated
#' @param AWX AWX dataframe
#' @param params list of parameters
#'
#' @return a function that computes likelihood for the the parameter which_not_known.
#' The function's input is always a vector where the parameters are ordered alphabetically by their names, gamma before lambda_0 before theta.
#' @export
twoKnownLik <- function(which_not_known,params, AWX){
  gamma <- params$gamma
  lambda_0 <- params$lambda_0
  theta <- params$theta

  negLik <- Vectorize(
    FUN = negLogLik,
    vectorize.args = which_not_known,
    SIMPLIFY = TRUE
  )

  if (which_not_known == "gamma") {
  oneLik <- purrr::partial(negLik,AWX = AWX, lambda_0 = lambda_0, theta = theta)
  }
  if (which_not_known == "lambda_0") {
    oneLik <- purrr::partial(negLik,AWX = AWX, gamma = gamma, theta = theta)

  }
  if (which_not_known == "theta") {
    oneLik <- purrr::partial(negLik,AWX = AWX, gamma = gamma ,lambda_0 = lambda_0)

  }
  return(Vectorize(oneLik))
}


#' Maximum likelihood estimation
#'
#' @param AWX AWX dataframe
#' @param params vector c(gamma, lambda_0, theta)
#' @importFrom stats optim
#' @return A list with elements: ans - the mle's, boundary = logical indicating whether solution is on boundary
#' @export
#'
mleFull <- function(AWX, params) {
  g_l0_t <- c(params$gamma,params$lambda_0,params$theta)
  opt <-
    optim(
      g_l0_t,
      # note that PARAMS is temporary
      fn = negLogLik.vec,
      lower = g_l0_t / 10,
      upper = g_l0_t * 10,
      method = "L-BFGS-B",
      gr = gradLogLik,
      AWX = AWX
    )
  is_boundary <- any((opt$par == g_l0_t / 10) |
                       (opt$par == g_l0_t * 10))

  ans <- opt$par
  names(ans) <- c("gamma", "lambda_0", "theta")
  return(list(ans = ans, boundary = is_boundary))
}


#' Liron MLE
#'
#' @param AWX data
#' @param acc accuracy for binary search
#'
#' @return A vector of theta and lambda estimates
#' @export
mleLiron <- function(AWX, acc = 1e-4) {
  #Lambda MLE given an estimator for theta (exponential patience)
  lambda.MLE <- function(theta.hat, A, W, X) {
    n <- length(W)
    a <-
      exp(-theta.hat * W[2:n]) - exp(-theta.hat * (W[1:(n - 1)] + X[1:(n - 1)]))
    b <-
      theta.hat * pmax(rep(0, n - 1), A[2:n] - W[1:(n - 1)] - X[1:(n -
                                                                     1)])
    lambda.hat <- n * theta.hat / sum(a + b)
    return(lambda.hat)
  }
  A <- AWX$A
  W <- AWX$W
  X <- AWX$X
  n <- length(W)
  Theta <- c(0, 10 ^ 3) #Search range
  d <- acc * 2
  #Bisection search for optimal theta
  while (abs(d) > acc)
  {
    theta.hat <- mean(Theta)
    lambda.hat <-
      lambda.MLE(theta.hat, A, W, X) #Lambda mle for given theta
    a <-
      (1 + theta.hat * W[2:n]) * exp(-theta.hat * W[2:n]) - (1 + theta.hat *
                                                               (W[1:(n - 1)] + X[1:(n - 1)])) * exp(-theta.hat * (W[1:(n - 1)] + X[1:(n -
                                                                                                                                        1)]))
    d <- mean(W[2:n]) - lambda.hat * mean(a) / (theta.hat ^ 2)
    #Update value:
    if (d > acc) {
      Theta[2] <- theta.hat
    }
    if (d < -acc) {
      Theta[1] <- theta.hat
    }
  }
  return(c("theta.liron" = theta.hat, "lambda.liron" = lambda.hat))
}


#' Mle with one known parameter
#'
#' @param AWX data
#' @param params parameters
#' @param which_known name of the _known_ parameter
#'
#' @return Maximum likelihood estimates when one known parameter
#' @export
mleOneKnown <- function(AWX, params, which_known){

  which_uknown <- setdiff(names(GL0T(params)), which_known)

  if (length(which_known) == 1){ # one uknown parameter
    knownLik <- oneKnownLik(which_known = which_known,params = params,AWX = AWX)
    # gradient to be implemented later
    opt <- optim(par = GL0T(params[which_uknown]),
          fn = knownLik)

    # KG/L/T for known gamma,lambda_0, theta, resepctively
    estimator_extension <- paste0("_K",toupper(which_known %>% stringr::str_sub(1,1)))
    mle <- opt$par
    names(mle) <- paste0(which_uknown, estimator_extension)
  }
  return(mle)
}

