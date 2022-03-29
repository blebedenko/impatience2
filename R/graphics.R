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


#' Plot the queue dynamics in a specified time interval
#'
#' @param from_to vector of start and stop of the time interval. Defaults to c(1,2).
#' @param RES full simulations results list
#' @param params scenario parameters
#'
#' @return plot only
#' @export
#'
#' @examples
#' R <- exampleDataRES()
#' pltQueueInterval(RES = R, params = exampleParams())
pltQueueInterval <- function(from_to = c(1,2), RES, params){


  A <- RES$Aj # start from 0
  Q <- RES$Qj
  A_tilde <- cumsum(A) # arrival times

  start_time <- from_to[1]
  end_time <- from_to[2]
  interval_indexes <- which(A_tilde >= start_time & A_tilde <= end_time)

  # subset only the intervals' obs
  A_tilde_interval <- A_tilde[interval_indexes]
  Q_interval <- Q[interval_indexes]
  # write into a dataframe for use with latticeExtra
  q_data <- data.frame(arrival_time = A_tilde_interval, queue = Q_interval)
  # obtain the rate function:
  rate <- rateFactory(params = params)
  # compute the rate for a grid over the interval
  time <- seq(start_time,end_time,length.out = 1000)
  r_data <- data.frame(time = time, rate = rate(time))

  p1 <-
    lattice::xyplot(
      queue ~ arrival_time,
      type = "s",
      lwd = 2,
      data = q_data,
      xlab = "time",
      col = "blue"
    )

  p2 <-
    lattice::xyplot(
      rate ~ time,
      data = r_data,
      type = "l",
      col = "purple",
      lwd = 3,
      lty = 3,
      key = list(
        #corner = c(0,0),
        space = c("bottom"),
        lines =
          list(
            col = c("purple", "blue"),
            lty = c(3, 1),
            lwd = c(3, 2)
          ),
        text =
          list(c(expression(lambda(t)),expression(Q(t))))
      )
    )
  latticeExtra::doubleYScale(p1, p2,add.ylab2 = T)
}



#' Plot the queue dynamics in a specified time interval
#'
#' @param segments number of segments to divide one cycle into. Defaults to 24 (hours per day).
#' @param RES full simulations results list
#' @param params scenario parameters
#' @export
pltArrivalsByTimeSegment <- function(segments = 24, RES,params){
  if ( length(segments) != 1 || segments <= 0)
    stop ("number of segments has to be a positive integer")
  A_tilde <- cumsum(RES$Aj)
  time_of_day <- A_tilde - trunc(A_tilde)
  breaks <- seq(0,1,by = 1 / segments)
  arrivals_segmented <- cut(time_of_day,breaks)
  levels(arrivals_segmented) <- 1:segments
  dat <- data.frame(time = time_of_day, queue = RES$Qj, waiting = pmax(RES$Qj - params$s, 0) )
  dat$segment <- cut(dat$time,breaks,labels = 1:segments)
  dat <- dat %>% group_by(segment) %>% summarise(L = mean(queue),
                                                 Lw = mean(waiting),
                                                 arrivals = n())
  p1 <- lattice::xyplot(arrivals ~ segment, data = dat, type = "h", lwd = 3, lty = 3, col = 1)
  p2 <- lattice::xyplot(dat$Lw~ as.numeric(segment) - 0.5, data = dat, type = "h", lty = 1, ylab = "No. Waiting", col = 2, lwd = 4,
                        key = list(lines = list(col = 1:2 , lwd = 3:4),
                                   text = list(lab = c('Arrivals/segment','No. waiting / segment')),
                                    title = "Arrivals and queue length"))
  latticeExtra::doubleYScale(p1,p2,add.ylab2 = T,under = )
}

