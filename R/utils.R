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
