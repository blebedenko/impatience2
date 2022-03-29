rm(list = ls())
library(devtools)
#imports <- c("base", "stats", "graphics", "utils", "svMisc", "svDialogs", "withr", "dplyr", "doParallel", "magrittr", "foreach")

RES <- exampleDataRES("C4")
params <- scenarioParams("C4")
pltQueueInterval(RES = RES, params = exampleParams(),from_to = c(5,10))
pltArrivalsByTimeSegment(RES = RES,params = params,segments = 24)
RES$Q
#' Plot the queue dynamics in a specified time interval
#'
#' @param segments number of segments to divide one cycle into. Defaults to 24 (hours per day).
#' @param RES full simulations results list
#' @param params scenario parameters
#' @export
pltArrivalsByTimeSegment <- function(segments = 24, RES,params){
  if (!is.integer(segments) || length(segments) != 1 || segments <= 0)
    stop ("number of segments has to be a positive integer")
  A_tilde <- cumsum(RES$Aj)
  time_of_day <- A_tilde - trunc(A_tilde)
  breaks <- seq(0,1,by = 1 / segments)
  levels(arrivals_segmented) <- 1:segments
  dat <- data.frame(time = time_of_day, queue = RES$Qj)
  dat$segment <- cut(dat$time,breaks,labels = 1:segments)
  dat <- dat %>% group_by(segment) %>% summarise(Q = mean(queue), arrivals = n())
  p1 <- lattice::xyplot(arrivals ~ segment, data = dat, type = "h", lwd = 3, lty = 3)
  p2 <- lattice::xyplot(Q ~ segment, data = dat, type = "h", lty = 1)
  latticeExtra::doubleYScale(p1,p2)
}

