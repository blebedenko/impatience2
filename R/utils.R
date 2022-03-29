#' Auxiliary function to create example parameters
#'
#' @return a list of example parameters
#' @export
#'
#' @examples
#' params <- listParams(gamma=1,lambda_0=2,theta=1, eta = 1, mu = 1 , s = 3)
#' identical(params,exampleParams())

exampleParams <- function() {
  ex_pars <-
    listParams(
      gamma = 1,
      lambda_0 = 2,
      theta = 1,
      eta = 1,
      mu = 1 ,
      s = 3
    )
  return(ex_pars)

}


#' Auxiliary function to create example dataset
#'
#' @return AWX dataset with n = 1000 generated with exampleParams()
#' @export
#'
#' @examples
#' set.seed(123)
#' dat <- exampleDataAWX()
#' set.seed(123)
#' pars <- exampleParams()
#' dat2 <- resSimAWX(n_thousands = 1,params = pars)
#' identical(dat,dat2)
exampleDataAWX <- function(){
  params <- exampleParams()
  AWX <- resSimAWX(n_thousands = 1,params = params)
  return(AWX)
}


#' get gamma, lambda_0 and theta
#'
#' @param params  names list of parameters (output of 'listParams()').
#'
#' @return a named vector (gamma,lambda_0,theta)
#' @export
#'
#' @examples
#' GL0T(exampleParams())
GL0T <-
  function(params){
    c(
      gamma = params$gamma,
      lambda_0 = params$lambda_0,
      theta = params$theta
    )
  }

#' Utility: turn RES to AWX
#'
#' @param RES results of resSimCosine()
#'
#' @return A data frame with columns A - inter-arrival times, W - waiting times, X - waiting time jump sizes
#' @export
#'
RES2AWX <-
  function(RES) {
    return(data.frame(
      A = RES$Aj,
      W = RES$Wj,
      X = RES$Xj
    ))
  }

#' parameters by simulation scenario
#'
#' @param scenario_name  "C1", "C2", "C3", "C4".
#'
#' @return params list with NA n. servers
#' @export
#'
#' @examples
#' scenarioParams("C3")
scenarioParams <- function(scenario_name){
  params <- switch (scenario_name,
    "C1" = listParams(gamma = 10,lambda_0 = 10,theta = 2,5,eta = 1,mu = 1, s=NA),
    "C2" = listParams(gamma = 40,lambda_0 = 10,theta = 2.5,eta = 1,mu = 1, s=NA),
    "C3" = listParams(gamma = 1,lambda_0 = 12,theta = 1,eta = 1,mu = 1, s=NA),
    "C4" = listParams(gamma = 100,lambda_0 = 50,theta = 10,eta = 1,mu = 1, s=NA)
  )
  return(params)
}


#' Pipe
#'
#' Put description here
#'
#' @importFrom magrittr %>%
#' @name %>%
#' @rdname pipe
#' @export
#' @param lhs,rhs specify what lhs and rhs are
#' @examples
#' # some examples if you want to highlight the usage in the package
NULL
