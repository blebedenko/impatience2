% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/simulation.R
\name{listParams}
\alias{listParams}
\title{Make Parameter list}
\usage{
listParams(gamma, lambda_0, theta, eta, mu, s, model = "cosine_exp")
}
\arguments{
\item{gamma}{periodical component of the rate function.}

\item{lambda_0}{constant component of the rate function.}

\item{theta}{exponential patience rate parameter.}

\item{eta}{shape parameter of the job size.}

\item{mu}{rate parameter of the job size.}

\item{s}{number of servers}

\item{model}{(default to "cosine_exp") names of the model to use}
}
\value{
named list with all the relevant parameters.
}
\description{
Make Parameter list
}
\examples{
listParams(gamma=10,lambda_0=20,theta=2.5, eta = 1, mu = 1 , s = 3)
}
