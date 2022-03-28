#' Plot the rate function
#'
#' @param params  names list of parameters (output of 'listParams()').
#' @param n_cycles of cycles to plot (defaults to 5)
#'
#' @export
#'
#' @examples
#' pltRate(exampleParams())
pltRate <- function(params, n_cycles = 5) {
  gamma <- params$gamma
  lambda_0 <- params$lambda_0
  lam <- function(t) {
    lambda_0 + (gamma / 2) * (cos(2 * pi * t) + 1)
  }
  graphics::curve(
    lam,
    from = 0,
    to = n_cycles,
    lwd = 2,
    xlab = "time",
    ylab = expression(lambda(t)),
    main =  bquote("Parameters:" ~ gamma == .(gamma) ~ "and" ~ lambda[0] == .(lambda_0))
  )
}
