#' (negative) mean log-likelihood for the sinusoidal arrivals model
#' The mean is returned instead of the sum - should help gradient based optimizers
#' @param gamma periodical component of rate function
#' @param lambda_0 constant component of rate function
#' @param theta parameter of exponential patience
#' @param AWX compact simulation results for memory economy (A,W,X).
#' @return The negative log-likelihood at the point provided. Based on mean (instead of sum).
#' @export
negLogLik <- function(gamma, lambda_0, theta, AWX) {
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
#' @return negative mean likelihood for the point g_l0_t
negLogLik.vec <- function(g_l0_t, AWX) {
  gamma <- g_l0_t[1]
  lambda_0 <- g_l0_t[2]
  theta <- g_l0_t[3]
  return(negLogLik(
      gamma = gamma,
      lambda_0 = lambda_0 ,
      theta = theta,
      AWX = AWX
    ))
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
gradLogLik <- function(gamma, lambda_0, theta, AWX) {
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

  names(negativeGradientMean) <- c("gamma", "lambda_0", "theta")

  return(negativeGradientMean)




}

#' Gradient of  likelihood, vector form
#'
#' @param g_l0_t named vector c(gamma, lambda_0, theta)
#' @param AWX AWX
#'
#' @return the gradient at (gamma,lambda_0,theta) for data AWX
#' @export
gradLogLik.vec <- function(g_l0_t, AWX) {
  gamma <- g_l0_t["gamma"]
  lambda_0 <- g_l0_t["lambda_0"]
  theta <- g_l0_t["theta"]
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

  names(negativeGradientMean) <- c("gamma", "lambda_0", "theta")

  return(negativeGradientMean)




}

#' generic gradient function
#'
#' @param point named vector of 2-3 of parameters (gamma,lambda_0,theta)
#' @param AWX the dataset
#' @param params the actual parameter values
#'
#' @return a named vector with the value of the gradient at point, assuming the
#' unspecified parameter is \emph{known!}
#' @export
#'
#' @examples
#'
#' params <- exampleParams()
#' names(GL0T(params))
#' AWX <- exampleDataAWX()
#' point <- GL0T(params)[c("lambda_0","theta")]
#' grad(point = point,AWX = AWX,params = params)
#' point <- GL0T(params)[c("lambda_0","gamma")]
#' grad(point = point,AWX = AWX,params = params)
#' point <- GL0T(params)[c("gamma","lambda_0")]
#' grad(point = point,AWX = AWX,params = params)
#' point <- GL0T(params) + runif(3)
#' grad(point = point,AWX = AWX,params = params)
grad <- function(point,params,AWX){
  if( length(point) == 3){ # no parameter is known
    return(gradLogLik.vec(g_l0_t = point, AWX = AWX))
  } else{
  # which_known <- setdiff(names(GL0T(params)),names(point)) # name of the known parameter
    return(gradOneKnown(point = point,params = params,AWX = AWX))
  }
}

#' gradient when one parameter is known
#'
#' @param point named vector of the \emph{uknown} parameters.
#' @param params simulation parameters
#' @param AWX dataset
#'
#' @return The gradient at the supplied point
#' @export
#'
gradOneKnown <- function(point,params, AWX){
  which_known <- setdiff(names(GL0T(params)),names(point)) # name of the known parameter
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

  if (which_known == "gamma") {
    # set gamma to the true value:
    gamma <- params$gamma
    lambda_0 <- point["lambda_0"]
    theta <- point["theta"]
    # compute derivatives by lambda_0 and theta only
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


      # compute the negative of the gradient elements' mean
      negativeGradientMean <-
        -c(mean(dl_lambda_0), mean(dl_theta))

      names(negativeGradientMean) <- c("lambda_0", "theta")

    }
  if (which_known == "lambda_0") {
    theta <- point["theta"]
    gamma <- point["gamma"]
    # set lambda_0 to the true value:
    lambda_0 <- params$lambda_0
    # compute derivatives by gamma and theta only
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


    # compute the negative of the gradient elements' mean
    negativeGradientMean <-
      -c(mean(dl_gamma), mean(dl_theta))

    names(negativeGradientMean) <- c("gamma", "theta")

  }
  if (which_known == "theta") {
    lambda_0 <- point["lambda_0"]
    gamma <- point["gamma"]
    # set theta to the true value:
    theta <- params$theta
    # compute derivatives by gamma and theta only
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



      # compute the negative of the gradient elements' mean
      negativeGradientMean <-
        -c(mean(dl_gamma), mean(dl_lambda_0))

      names(negativeGradientMean) <- c("gamma", "lambda_0")
  }
  return(negativeGradientMean)


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
oneKnownLik <- function(which_known, params, AWX) {


  if (which_known == "gamma") {
    # negLik <- Vectorize(
    #   FUN = negLogLik,
    #   vectorize.args = c("lambda_0", "theta"),
    #   SIMPLIFY = TRUE
    # )
    # likelihood for lambda_0, theta with known gamma
    negLik <- purrr::partial(negLogLik, AWX = AWX, gamma = params$gamma)
    oneLik <- function(lambda_0, theta)
      negLik(lambda_0, theta)

  }
  if (which_known == "lambda_0") {
    negLik <- Vectorize(
      FUN = negLogLik,
      vectorize.args = c("gamma", "theta"),
      SIMPLIFY = TRUE
    )
    # likelihood for gamma, theta with known lambda_0
    negLik <- purrr::partial(negLik, AWX = AWX, lambda_0 = params$lambda_0)
    oneLik <- function(gamma, theta)
      negLik(gamma, theta)

  }
  if (which_known == "theta") {
    negLik <- Vectorize(
      FUN = negLogLik,
      vectorize.args = c("gamma", "lambda_0"),
      SIMPLIFY = TRUE
    )
    # likelihood for gamma, lambda_0 with known theta
    negLik <- purrr::partial(negLik, AWX = AWX, theta = params$theta)
    oneLik <- function(gamma, lambda_0)
      negLik(gamma, lambda_0)
  }

  curr_params <- setdiff(c("gamma", "lambda_0", "theta"), which_known)

  negLogL <-
    function(two_par_vec)
      oneLik(two_par_vec[1], two_par_vec[2]) # wrap it for optim
  return(negLogL)
}

#' Title f
#'
#' @inheritParams oneKnownLik
#' @return a gradient function for the two other parameter
#' @export
oneKnownGrad <- function(which_known, params, AWX) {
  # gamma <- params$gamma
  # lambda_0 <- params$lambda_0
  # theta <- params$theta
  which_uknown <- setdiff(names(GL0T(params)), which_known)
  if (which_known == "gamma") {
    # gradient for lambda_0, theta with known gamma
    oneGrad <-
      purrr::partial(gradLogLik, AWX = AWX, gamma = params$gamma)

  }

  if (which_known == "lambda_0") {
    oneGrad <-
      purrr::partial(gradLogLik, AWX = AWX, lambda_0 = params$lambda_0)

  }
  if (which_known == "theta") {
    oneGrad <-
      purrr::partial(gradLogLik, AWX = AWX, theta = params$theta)

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
twoKnownLik <- function(which_not_known, params, AWX) {
  gamma <- params$gamma
  lambda_0 <- params$lambda_0
  theta <- params$theta

  negLik <- Vectorize(FUN = negLogLik,
                      vectorize.args = which_not_known,
                      SIMPLIFY = TRUE)

  if (which_not_known == "gamma") {
    oneLik <-
      purrr::partial(negLogLik,
                     AWX = AWX,
                     lambda_0 = lambda_0,
                     theta = theta)
  }
  if (which_not_known == "lambda_0") {
    oneLik <-
      purrr::partial(negLogLik,
                     AWX = AWX,
                     gamma = gamma,
                     theta = theta)

  }
  if (which_not_known == "theta") {
    oneLik <-
      purrr::partial(negLogLik,
                     AWX = AWX,
                     gamma = gamma ,
                     lambda_0 = lambda_0)

  }
  return(Vectorize(oneLik))
}



#' Compute all MLE's for a dataset
#'
#' @param AWX can be a dataframe or a path
#' @param params parameter list
#'
#' @return named vector of all MLE's
#' @export
#'
mleAll <- function(AWX, params) {
  # if data provided
  if (is.data.frame(AWX)) {
    mle <- mleGetAll(AWX = AWX, params = params)
  }
  if (is.character(AWX)) {
      mle <- mleGetAll(AWX = utils::read.csv(AWX), params = params)

  }
  return(mle)
}

#' get all MLE's for dataset
#'
#' @param AWX dataframe
#' @param params parameter list
#'
#' @return a named vector with all possible MLE's
#' @export
#'
mleGetAll <- function(AWX, params) {
  mle <- c(
    mleFull(AWX, params),
    mleLiron(AWX),
    mleOneKnown(
      AWX = AWX,
      params = params,
      which_known = "gamma"
    ),
    mleOneKnown(
      AWX = AWX,
      params = params,
      which_known = "lambda_0"
    ),
    mleOneKnown(
      AWX = AWX,
      params = params,
      which_known = "theta"
    ),
    mleTwoKnown(
      AWX = AWX,
      params = params,
      which_not_known = "gamma"
    ),
    mleTwoKnown(
      AWX = AWX,
      params = params,
      which_not_known = "lambda_0"
    ),
    mleTwoKnown(
      AWX = AWX,
      params = params,
      which_not_known = "theta"
    )
  )

  return(mle)
}


#' Maximum likelihood estimation - Full model
#'
#' @param AWX AWX dataframe
#' @param params vector c(gamma, lambda_0, theta)
#' @importFrom stats optim
#' @return A list with elements: ans - the mle's, boundary = logical indicating whether solution is on boundary
#' @export
#'
mleFull <- function(AWX, params) {
  g_l0_t <- GL0T(params)
  opt <-
    optim(
      g_l0_t,
      # note that PARAMS is temporary
      fn = negLogLik.vec,
      lower = g_l0_t / 10,
      upper = g_l0_t * 10,
      method = "L-BFGS-B",
      gr = gradLogLik.vec,
      AWX = AWX
    )
  is_boundary <- any((opt$par == g_l0_t / 10) |
                       (opt$par == g_l0_t * 10))

  ans <- opt$par
  names(ans) <- c("gamma", "lambda_0", "theta")
  return(c(boris = ans, boundary = is_boundary))
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
mleOneKnown <- function(AWX, params, which_known) {
  which_uknown <- setdiff(names(GL0T(params)), which_known)

  if (length(which_known) == 1) {
    # one uknown parameter
    knownLik <-
      oneKnownLik(which_known = which_known,
                  params = params,
                  AWX = AWX)
    # gradient to be implemented later
    opt <- optim(par = GL0T(params[which_uknown]),
                 upper = GL0T(params[which_uknown]) * 10,
                 lower = GL0T(params[which_uknown]) / 10,
                 method = "L-BFGS-B",
                 fn = knownLik)

    # KG/L/T for known gamma,lambda_0, theta, resepctively
    estimator_extension <-
      paste0("_K", toupper(which_known %>% stringr::str_sub(1, 1)))
    mle <- opt$par
    names(mle) <- paste0(which_uknown, estimator_extension)
  }
  return(mle)
}



#' Mle with two known parameters
#'
#' @param AWX data
#' @param params parameters
#' @param which_not_known name of the UNKNOWN parameter to be estimated
#'
#' @return Maximum likelihood estimates when one known parameter
#' @export
mleTwoKnown <- function(AWX, params, which_not_known) {
  gamma <- params$gamma
  lambda_0 <- params$lambda_0
  theta <- params$theta

  which_known <- setdiff(names(GL0T(params)), which_not_known)
  if (which_not_known == "gamma") {
    oneLik <-
      purrr::partial(negLogLik,
                     AWX = AWX,
                     lambda_0 = lambda_0,
                     theta = theta)
  }
  if (which_not_known == "lambda_0") {
    oneLik <-
      purrr::partial(negLogLik,
                     AWX = AWX,
                     gamma = gamma,
                     theta = theta)

  }
  if (which_not_known == "theta") {
    oneLik <-
      purrr::partial(negLogLik,
                     AWX = AWX,
                     gamma = gamma ,
                     lambda_0 = lambda_0)

  }
  oneLik.vec <- Vectorize(oneLik)
  opt <- stats::optimize(f = oneLik.vec,
                         lower = params[[which_not_known]] / 10,
                         upper = params[[which_not_known]] * 10)
  mle <- opt$minimum
  names(mle) <- paste0(which_not_known, "_only")
  return(mle)
}
