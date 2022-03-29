#'  Make Parameter list
#'
#' @param gamma periodical component of the rate function.
#' @param lambda_0  constant component of the rate function.
#' @param theta  exponential patience rate parameter.
#' @param eta shape parameter of the job size.
#' @param mu rate parameter of the job size.
#' @param s number of servers
#' @param model (default to "cosine_exp") names of the model to use
#' @return named list with all the relevant parameters.
#' @export
#'
#' @examples
#' listParams(gamma=10,lambda_0=20,theta=2.5, eta = 1, mu = 1 , s = 3)
listParams <-
  function(gamma, lambda_0, theta, eta, mu, s,  model = "cosine_exp") {
    par_list <-
      list(
        gamma = gamma,
        lambda_0 = lambda_0,
        theta = theta,
        eta = eta,
        mu = mu,
        s = s
      )
    return(par_list)
  }



#' Rate function factory
#'
#' @param params  named list of parameters (output of 'listParams()').
#'
#' @return A function with argument t that computes the arrival rate.
#' @export
#'
#' @examples
#' params <- listParams(gamma=10,lambda_0=20,theta=2.5, eta = 1, mu = 1 , s = 3)
#' rate <- rateFactory(params)
#' rate(t = runif(3))
rateFactory <- function(params) {
  gamma <- params$gamma
  lambda_0 <- params$lambda_0
  return(function(t)
    lambda_0 + (gamma / 2) * (1 + cos(2 * pi * t)))
}



#' Next Arrival of inhomogeneous Poisson process
#'
#' Generates  next arrival time (not inter-arrival!) of a nonhomogeneous Poisson process.
#' @param current_time the time on the clock.
#' @param params  named list of parameters (output of 'listParams()').

#' @return time of the next arrival.
#'
#' @export
#'
#' @examples
#' params <- listParams(gamma=10,lambda_0=20,theta=2.5, eta = 1, mu = 1 , s = 3)
#' nextArrivalCosine(1.2,params = params)
nextArrivalCosine <- function(current_time, params) {
  gamma <- params$gamma
  lambda_0 <- params$lambda_0
  lambda_sup <- gamma + lambda_0 # highest rate in the cosine model
  rate <- rateFactory(params)
  arrived <- FALSE
  while (!arrived) {
    u1 <- stats::runif(1)
    current_time <-
      current_time - (1 / lambda_sup) * log(u1) # generate next arrival from Pois(sup_lambda)
    u2 <- stats::runif(1)
    arrived <- u2 <= rate(current_time) / lambda_sup
  }
  return(current_time)
}

#' Customer's job size and patience
#' Provides with a realization of Y~Exp(theta) and B~Gamma(eta,mu)
#' @param params  named list of parameters (output of 'listParams()').
#'
#' @return list with elements 'Patience' and 'Jobsize'
#' @export
#' @examples
#' params <- listParams(gamma=10,lambda_0=20,theta=2.5, eta = 1, mu = 1 , s = 3)
#'  customer <- customerExp(params)
#'  customer
customerExp <- function(params) {
  theta <- params$theta
  eta <- params$eta
  mu <- params$mu
  B <-  stats::rgamma(1, shape = eta, rate = mu) #Job sizes
  Y <- stats::rexp(1, theta)
  return(list(Patience = Y, Jobsize = B))
}


#' Liron's Virtual Waiting function
#'
#' @param Res.service vector of residual service times
#' @param s the number of servers
#' @export
#' @examples
#' VW(c(1.2,1.5,0.3), s = 2)
#' VW(c(1.2,1.5,0.3), s = 3)
#' VW(c(1.2,1.5,0.3), s = 4)
VW <- function(Res.service, s) {
  Q.length <- length(Res.service) #Number in the system
  if (Q.length < s) {
    virtual.wait <- 0
  } else
  {
    D.times <- rep(NA, Q.length + 1)
    Vw <- rep(NA, Q.length + 1)
    D.times[1:s] <- Res.service[1:s]
    Vw[1:s] <- rep(0, s) #VW for customers in service
    for (i in (s + 1):(Q.length + 1))
    {
      D.i <- sort(D.times[1:i]) #Sorted departures of customers ahead of i
      Vw[i] <- D.i[i - s] #Virtual waiting time for position i
      if (i <= Q.length) {
        D.times[i] <- Res.service[i] + Vw[i]
      } #Departure times
    }
    virtual.wait <- Vw[Q.length + 1]
  }
  return(virtual.wait)
}



#' Atomic Simulation for periodic arrivals (cosine model)
#' @param n Number of samples to generate.
#' @param params  named list of parameters (output of 'listParams()').
#' @return A list with the queue data
#' @export
#' @examples
#' params <- exampleParams()
#' R <- resSimCosine(n=100, params = params)
resSimCosine <- function(n, params) {
  gamma <- params$gamma
  lambda_0 <- params$lambda_0
  theta <- params$theta
  eta <- params$eta
  mu <- params$mu
  s <- params$s

  #find the supremum of the arrival function
  lambda_sup <- lambda_0 + gamma
  lambdaFunction <-   rateFactory(params)
  #Running simulation variables
  klok <- 0
  n.counter <- 0 #Observation counter
  Q.length <- 0 #Queue length (including service)
  virtual.wait <- 0 #Virtual waiting time process
  m <- 0 #Total customer counter
  #A <- rexp(1,lambda) #First arrival time
  A <- nextArrivalCosine(klok[length(klok)], params = params)
  A <- A - klok[length(klok)]
  klok <- c(klok, klok[length(klok)] + A)
  trans.last <- A #Counter of time since last transition
  time.last <- A #Counter of time since last admission event
  Res.service <- numeric(0) #Vector of residual service times

  #Output vectors:
  Wj <- rep(NA, n) #Vector of workloads before jumps
  Xj <- rep(NA, n) #Vector of workload jumps
  Aj <- rep(NA, n) #Vector of effective inter-arrival times
  Qj <- rep(NA, n) #Vector of queue lengths at arrival times
  IPj <- A #Vector of idle periods
  Yj <-
    rep(NA, n) # Vector of patience values of customers that join
  Q.trans <- 0 #Vector of queue lengths at transitions
  IT.times <- numeric(0) #Vector of inter-transition times
  Pl <- 1 #Proportion of lost customers
  Nl <- 0 # number of lost customers

  while (n.counter < n + 1)
  {
    m <- m + 1
    customer <- customerExp(params = params)
    #Generate patience and job size:
    B <- customer$Jobsize
    Y <- customer$Patience
    if (virtual.wait <= Y)
      #New customer if patience is high enough
    {
      n.counter <- n.counter + 1 #Count observation
      Res.service <-
        c(Res.service, B) #Add job to residual service times vector
      Q.length <- length(Res.service) #Queue length
      Q.trans <- c(Q.trans, Q.length) #Add current queue length
      #Add new observations:
      Wj[n.counter] <- virtual.wait
      Aj[n.counter] <- time.last
      Qj[n.counter] <-
        Q.length - 1 #Queue length (excluding new arrival)
      Yj[n.counter] <- Y # patience of the customer arriving
      IT.times <- c(IT.times, trans.last) #Update transition time
      trans.last <- 0 #Reset transition time
      time.last <- 0 #Reset last arrival time
      Pl <- m * Pl / (m + 1) #Update loss proportion
    } else {
      Pl <- m * Pl / (m + 1) + 1 / (m + 1)
      Nl <- Nl + 1
    }

    #Update system until next arrival event
    #A <- rexp(1,lambda) #Next arrival time
    A <- nextArrivalCosine(klok[length(klok)], params = params)
    A <- A - klok[length(klok)]
    klok <- c(klok, klok[length(klok)] + A)
    time.last <-
      time.last + A #Add arrival time to effective arrival time

    #Departure and residual service times of customers in the system:
    Q.length <- length(Res.service) #Queue length
    D.times <- rep(NA, Q.length) #Departure times
    Vw <- rep(NA, Q.length) #Virtual waiting times
    for (i in 1:Q.length)
    {
      if (i <= s)
      {
        Vw[i] <- 0 #No virtual waiting time
        D.times[i] <-
          Res.service[i] #Departure time is the residual service time
        Res.service[i] <-
          max(Res.service[i] - A, 0) #Update residual service time
      } else
      {
        D.i <- sort(D.times[1:i]) #Sorted departures of customers ahead of i
        Vw[i] <- D.i[i - s] #Time of service start for customer i
        D.times[i] <- Res.service[i] + Vw[i] #Departure time
        serv.i <-
          max(0, A - Vw[i]) #Service obtained before next arrival
        Res.service[i] <-
          max(Res.service[i] - serv.i, 0) #New residual service
      }
    }
    #Jump of virtual waiting time:
    if (virtual.wait <= Y)
    {
      if (Q.length < s)
      {
        Xj[n.counter] <- 0
      } else
      {
        Xj[n.counter] <- sort(D.times)[Q.length + 1 - s] - virtual.wait
      }
    }
    #Update residual service times:
    Res.service <-
      Res.service[!(Res.service == 0)] #Remove completed services

    #Update transition times and queue lengths:
    D.before <-
      which(D.times <= A) #Departing customers before next arrival
    if (length(D.before) > 0)
    {
      T.d <- sort(D.times[D.before]) #Sorted departure times
      for (i in 1:length(D.before))
      {
        Q.trans <-
          c(Q.trans, Q.length - i) #Update queue length at departures
        if (i == 1)
        {
          trans.last <- trans.last + T.d[1] #Update time since last transition
          IT.times <-
            c(IT.times, trans.last) #Departure transition time
          trans.last <- 0 #Reset transition time
        } else
        {
          trans.last <-
            trans.last + T.d[i] - T.d[i - 1] #Update time since last transition
          IT.times <-
            c(IT.times, trans.last) #Departure transition time
          trans.last <- 0 #Reset transition time
        }
        if (Q.trans[length(Q.trans)] == 0) {
          IPj <- cbind(IPj, A - T.d[length(T.d)])
        } #Add idle time observation
      }
      trans.last <-
        A - T.d[i] #Update remaining time until next arrival
    } else if (length(D.before) == 0)
    {
      trans.last <-
        trans.last + A #Update timer since least transition with new arrival
    }
    virtual.wait <- VW(Res.service, s) #Update virtual waiting time

  }
  # progress bar
  if (m %% 1000 == 0)
    cat("* ")

  RES <- list(
    Aj = Aj,
    # interarrival times
    Xj = Xj,
    # worlkload jumps upon arrival
    Wj = Wj,
    # waiting time experienced
    Qj = Qj,
    # queue length at arrival
    IPj = IPj,
    # idle periods
    Q.trans = Q.trans,
    #queue @ transitions
    IT.times = IT.times,
    #intertransition times
    Yj = Yj,
    # patiences
    klok = klok,
    # the clock ticks
    Pl = Pl,
    # proportion lost customers
    Nl = Nl,
    # number of lost customers
    Res.service = Res.service,
    # residual service

    virtual.wait = virtual.wait # virtual waiting

  )

  return(RES)
}




#' Simulate results, possibly from given initial state
#'
#' @param E0 Optional argument - A list with results from a previous simulation
#' @param n Number of samples to generate.
#' @param params  named list of parameters (output of 'listParams()').
#'
#' @return same as resSimCosine
#' @export
#'
#' @examples
#'  params <- exampleParams()
#' R <- resSimCosine(n=100, params = params)
#' R2 <- resSimCosineContinue(E0 = R, n = 100, params = params)
resSimCosineContinue <-
  function(E0 = NULL, n, params) {
    # auxiliary function:
   # parameters:
     gamma <- params$gamma
    lambda_0 <- params$lambda_0
    theta <- params$theta
    eta <- params$eta
    mu <- params$mu
    s <- params$s
    #find the supremum of the arrival function
    lambda_sup <- lambda_0 + gamma
    lambdaFunction <-   rateFactory(params = params)
    s <- params$s
    # if no initial conditions provided, generate as usual:
    if (is.null(E0)) {
      RES <- resSimCosine(
        # the "regular" generation function
        n = n,
       params = params
      )
    } else {
      # meaning that initial conditions were provided

      # Running simulation variables - initialized from E0
      klok <- dplyr::last(E0$klok)
      n.counter <- 0 # Observation counter for these results
      Q.length <-
        length(E0$Res.service) # last known queue length (including service)
      virtual.wait <-
        E0$virtual.wait #Virtual waiting time process
      m <- 0 # Total customer counter
      A <- nextArrivalCosine(klok[length(klok)],  params = params)
      A <- A - klok[length(klok)]
      klok <- c(klok, klok[length(klok)] + A)
      trans.last <- A #Counter of time since last transition
      time.last <- A #Counter of time since last admission event
      Res.service <-
        E0$Res.service #Vector of residual service times

      #Output vectors:
      Wj <- rep(NA, n) #Vector of workloads before jumps
      Xj <- rep(NA, n) #Vector of workload jumps
      Aj <- rep(NA, n) #Vector of effective inter-arrival times
      Qj <- rep(NA, n) #Vector of queue lengths at arrival times
      IPj <- A #Vector of idle periods
      Yj <-
        rep(NA, n) # Vector of patience values of customers that join
      Q.trans <-
        dplyr::last(E0$Q.trans) #Vector of queue lengths at transitions
      IT.times <- numeric(0) #Vector of inter-transition times
      Pl <- 1 #Proportion of lost customers
      Nl <- 0 # number of lost customers

      while (n.counter < n + 1)
      {
        m <- m + 1
        customer <- customerExp(params = params)
        #Generate patience and job size:
        B <- customer$Jobsize
        Y <- customer$Patience
        if (virtual.wait <= Y)
          #New customer if patience is high enough
        {
          n.counter <- n.counter + 1 #Count observation
          Res.service <-
            c(Res.service, B) #Add job to residual service times vector
          Q.length <- length(Res.service) #Queue length
          Q.trans <-
            c(Q.trans, Q.length) #Add current queue length
          #Add new observations:
          Wj[n.counter] <- virtual.wait
          Aj[n.counter] <- time.last
          Qj[n.counter] <-
            Q.length - 1 #Queue length (excluding new arrival)
          Yj[n.counter] <- Y # patience of the customer arriving
          IT.times <-
            c(IT.times, trans.last) #Update transition time
          trans.last <- 0 #Reset transition time
          time.last <- 0 #Reset last arrival time
          Pl <- m * Pl / (m + 1) #Update loss proportion
        } else {
          Pl <- m * Pl / (m + 1) + 1 / (m + 1)
          Nl <- Nl + 1
        }

        #Update system until next arrival event
        #A <- rexp(1,lambda) #Next arrival time
        A <-
          nextArrivalCosine(klok[length(klok)],  params = params)
        A <- A - klok[length(klok)]
        klok <- c(klok, klok[length(klok)] + A)
        time.last <-
          time.last + A #Add arrival time to effective arrival time

        #Departure and residual service times of customers in the system:
        Q.length <- length(Res.service) #Queue length
        D.times <- rep(NA, Q.length) #Departure times
        Vw <- rep(NA, Q.length) #Virtual waiting times
        for (i in 1:Q.length)
        {
          if (i <= s)
          {
            Vw[i] <- 0 #No virtual waiting time
            D.times[i] <-
              Res.service[i] #Departure time is the residual service time
            Res.service[i] <-
              max(Res.service[i] - A, 0) #Update residual service time
          } else
          {
            D.i <- sort(D.times[1:i]) #Sorted departures of customers ahead of i
            Vw[i] <-
              D.i[i - s] #Time of service start for customer i
            D.times[i] <- Res.service[i] + Vw[i] #Departure time
            serv.i <-
              max(0, A - Vw[i]) #Service obtained before next arrival
            Res.service[i] <-
              max(Res.service[i] - serv.i, 0) #New residual service
          }
        }
        #Jump of virtual waiting time:
        if (virtual.wait <= Y)
        {
          if (Q.length < s)
          {
            Xj[n.counter] <- 0
          } else
          {
            Xj[n.counter] <- sort(D.times)[Q.length + 1 - s] - virtual.wait
          }
        }
        #Update residual service times:
        Res.service <-
          Res.service[!(Res.service == 0)] #Remove completed services

        #Update transition times and queue lengths:
        D.before <-
          which(D.times <= A) #Departing customers before next arrival
        if (length(D.before) > 0)
        {
          T.d <- sort(D.times[D.before]) #Sorted departure times
          for (i in 1:length(D.before))
          {
            Q.trans <-
              c(Q.trans, Q.length - i) #Update queue length at departures
            if (i == 1)
            {
              trans.last <- trans.last + T.d[1] #Update time since last transition
              IT.times <-
                c(IT.times, trans.last) #Departure transition time
              trans.last <- 0 #Reset transition time
            } else
            {
              trans.last <-
                trans.last + T.d[i] - T.d[i - 1] #Update time since last transition
              IT.times <-
                c(IT.times, trans.last) #Departure transition time
              trans.last <- 0 #Reset transition time
            }
            if (Q.trans[length(Q.trans)] == 0) {
              IPj <- cbind(IPj, A - T.d[length(T.d)])
            } #Add idle time observation
          }
          trans.last <-
            A - T.d[i] #Update remaining time until next arrival
        } else if (length(D.before) == 0)
        {
          trans.last <-
            trans.last + A #Update timer since least transition with new arrival
        }
        virtual.wait <-
          VW(Res.service, s) #Update virtual waiting time

      }
      # progress bar
      if (m %% 1000 == 0)
        cat("* ")



      RES <- list(
        Aj = Aj,
        # interarrival times
        Xj = Xj,
        # worlkload jumps upon arrival
        Wj = Wj,
        # waiting time experienced
        Qj = Qj,
        # queue length at arrival
        IPj = IPj,
        # idle periods
        Q.trans = Q.trans,
        #queue @ transitions
        IT.times = IT.times,
        #intertransition times
        Yj = Yj,
        # patiences
        klok = klok,
        # the clock ticks
        Pl = Pl,
        # proportion lost customers
        Nl = Nl,
        # number of lost customers
        Res.service = Res.service,
        # residual service
        virtual.wait = virtual.wait # virtual waiting

      )
    }

    return(RES)
  }





#' Simulation of big sample sizes, by parts
#'
#' @param n_thousands what is the desired n (in thousands)?
#' @param params  named list of parameters (output of 'listParams()').
#' @return AWX dataframe
#' @export
#'
#' @examples
#' params <-  exampleParams()
#' a <- Sys.time()
#' res <- resSimAWX(3, params = params)# 10 K
#' Sys.time() - a
resSimAWX <-
  function(n_thousands,
           params) {
    gamma <- params$gamma
    lambda_0 <- params$lambda_0
    theta <- params$theta
    eta <- params$eta
    mu <- params$mu
    s <- params$s
    n_obs <- 1000 # always generate by pieces of 1000 arrivals
    n <- n_obs # to keep consistency

    # The first results:
    initial_RES <- resSimCosine(
      n = n_obs,
     params = params
    )

    d1 <- RES2AWX(initial_RES)
    last_RES <- initial_RES
    # the other simulations:

    for (i in 2:n_thousands) {
      next_RES <- resSimCosineContinue(
        E0 = last_RES,
        n = n_obs,
       params = params
      )
      svMisc::progress(i, n_thousands)
      # if (i %% 10 == 0)

      d2 <- RES2AWX(next_RES) # take the data
      d1 <- rbind(d1, d2) # run over the d1 data
      last_RES <- next_RES
    }
    # as each iteration produces 1001 - we omit the last "extra" observations
    return(utils::head(d1,n_thousands * 1000L))

  }



#' Interactive function that creates an entire scenario
#'
#' @return
#' @details User is prompted for the parameters in the console.
#' @importFrom foreach `%do%`
#' @importFrom foreach `%dopar%`
#' @export
#' @return Returns NULL. Creates folders with the realizations for each value of s.
makeAWXDirectories <- function() {

  withr::local_dir(svDialogs::dlg_dir()$res) # set the directory to where pointed
  n_cores <-
    as.numeric(readline(prompt = "How many cores to use? "))
  n_thousands <-
    as.numeric(readline(prompt = "n (thousands) = : "))
  N_files <- as.numeric(readline(prompt = "How many files? "))
  lambda_0 <-
    as.numeric(readline(prompt = "Input a value for lambda_0: "))
  gamma <-
    as.numeric(readline(prompt = "Input a value for gamma: "))
  theta <-
    as.numeric(readline(prompt = "Input a value for theta: "))
  eta <- as.numeric(readline(prompt = "Input a value for eta: "))
  mu <- as.numeric(readline(prompt = "Input a value for mu: "))
  prompt <-
    "enter numbers of servers for the experiments (space-separated) \n"
  s_values <- as.integer(strsplit(readline(prompt), " ")[[1]])

  min_rate <- lambda_0
  max_rate <- lambda_0 + gamma # not divided by two!
  ave_rate <- lambda_0 + gamma / 2
  rates <-
    c("minimum" = min_rate,
      "average" = ave_rate,
      "maximum" = max_rate)
  offered_loads <-
    foreach::foreach(s = s_values, .combine = rbind) %do% {
      rhos <- c(rates["minimum"] * eta / (s * mu),
                rates["average"] * eta / (s * mu),
                rates["maximum"] * eta / (s * mu))

    }

  all_params <- c(lambda_0, gamma, theta, eta, mu)
  param_names <- c("lambda_0", "gamma", "theta", "eta", "mu")
  param_message <-
    paste(param_names, " = ", all_params, "\n", collapse = "")
  server_message <-
    paste0("no. of servers: ", paste0(s_values, collapse = ","))
  obs_message <- paste0("no. observations: ", n_thousands,"K", "\n",
                        "no. of files: ", N_files)
  user_message <-
    paste0(param_message, "\n", server_message, "\n", obs_message)
  user_confirm <-
    svDialogs::dlg_message(message = user_message, type = "okcancel")$res
  user_confirm <- user_confirm == "ok"

  # message about offered loads:
  offered_loads <- data.frame(offered_loads)

  cat("The offered loads at the parameters you suggest are:\n")
  rownames(offered_loads) <- paste0("s=", s_values)
  print(offered_loads)
  final_confirm <-
    svDialogs::dlg_message(message = "Are you sure? To start, hit OK",
                           type = "okcancel")$res
  final_confirm <- final_confirm == "ok"
  if (!final_confirm)
    stop("Try again")

  # The actual computation:
  if (final_confirm)
  {
    cl <- parallel::makeCluster(n_cores)
    doParallel::registerDoParallel(cl)
    for (s in s_values) {
      # name for the current dir
      curr_dirname <-
        paste0("realizations for s=", s, "/")
      # create the dir
      dir.create(path = curr_dirname)
      # list the parameters anew for each value of s
      params <- listParams(gamma=gamma,lambda_0=lambda_0,theta=theta,eta=eta,mu=mu,s=s)
      cat(paste0(as.character(Sys.time()), " Starting now.\n", collapse = " "))
      # execute in the created directory

      foreach::foreach(
        i = 1:N_files,
        .packages = c("tidyverse", "impatience2","withr"),
        .combine = list
      ) %dopar% {
        RES <- resSimAWX(n_thousands = n_thousands,
                         params = params)

        dat <-
          RES # resSimAWX returns only AWX currently
        name <- filenamer(n_obs = n_thousands * 1000L,
                          params = params)

        withr::with_dir(curr_dirname, utils::write.csv(dat, file = name, row.names = FALSE))

      }




      cat("done with s = ", s, "...", "\n")

      gc()
    }


  }
  parallel::stopCluster(cl)
  cat(paste0(as.character(Sys.time()), " Finished!\n", collapse = " "))

}

#' Utility: name the individual simulation results
#'
#' @param n_obs number of observations in data
#' @param params  named list of parameters (output of 'listParams()').
#'
#' @return a time-stamped filename for a simulation dataframe
#' @export
#'
#' @examples
#' filenamer(1000,params = exampleParams())
filenamer <- function(n_obs, params) {
  gamma <- params$gamma
  lambda_0 <- params$lambda_0
  theta <- params$theta
  eta <- params$eta
  mu <- params$mu
  s <- params$s
  name <- (
    paste0(
      "AWX_",
      "s=",
      s,
      "n=",
      n_obs,
      "gamma=",
      gamma,
      "lambda_0=",
      lambda_0,
      "theta=",
      theta,
      "eta=",
      eta,
      "mu=",
      mu
    )
  )
  name <- paste0(name,
                 as.numeric(Sys.time()),
                 ".csv") # name with timestamp and extension
  return(name)
}

