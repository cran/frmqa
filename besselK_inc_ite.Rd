\name{besselK_inc_ite}
\alias{besselK_inc_ite}

\title{Approximations of the Incomplete BesselK Function Using an 
Iterative Algorithm}

\description{
    Approximates incomplete modified Bessel function using the algorithm
    provided by Slevinsky and Safouhi (2010) and gives warnings when 
    loss of accuracy or an exception occurs.
}
\usage{
    besselK_inc_ite(x, y, lambda, traceIBF = TRUE, epsilon = 0.95, nmax = 120)
}
\arguments{
    \item{x}{Argument x > 0.}
    \item{y}{Argument y > 0.}
    \item{lambda}{Order of the BesselK function.}
    \item{traceIBF}{Logical. Tracks the approximation process if \code{TRUE}.}
    \item{epsilon}{Determines approximation accuracy which equals machine 
    accuracy raised to the power of \code{epsilon}.}
    \item{nmax}{Integer. Maximum number of iterations required.}
}
\details{ One of the integral representations of the incomplete Bessel 
  function is given by \deqn{K_{\lambda}(x, y) = \int_1^{\infty}\,
  e^{-xt -y/t}\,t^{-\lambda-1}\,dt,} \code{besselK_inc_ite} approximates 
  this integral and gives warnings when encountering computational 
  problem or exceptions. 
}
\references{
    \cite{Slevinsky, R. M., and Safouhi, H (2010) A recursive algorithm 
    for the \eqn{G} transformation and accurate computation of 
    incomplete Bessel functions. \emph{Appl. Numer. Math.}, 
    \bold{60} 1411--1417}

    \cite{Function \code{incompleteBesselK} in package DistributionUtils}

    \cite{Kahan, W (1981), Why do we need a foating-point arithmetic standard?
        \url{http://www.cs.berkeley.edu/ wkahan/ieee754status/why-ieee.pdf}}
}
\note{
    A C++ version of this function will be added in due course.
}

\author{
  Thanh T. Tran \email{frmqa.package@gmail.com}
}        

\examples{
## "Exact" evaluation are given by Maple 15
options(digits = 16)

## Gives accurate approximations
## x = 0.01, y = 4, lambda = 0, exact value = 2.225 310 761 266 4692
besselK_inc_ite(0.01, 4, 0, traceIBF = FALSE, epsilon = 0.95, nmax = 160)
## x = 1, y = 4, lambda = 2, exact value = 0.006 101 836 926 254 8540
besselK_inc_ite(1, 4, 2, traceIBF = FALSE, epsilon = 0.95, nmax = 160)

## NaN occurs
## x = 1, y = 1.5, lambda = 8, exact value = 0.010 515 920 838 551 3164
\dontrun{besselK_inc_ite(1, 1.5, 8, traceIBF = TRUE, epsilon = 0.95, 
  nmax = 180)}

##Severe loss of accuracy
## x = 14.5, y = 19, lambda = 0, exact value = 9.556 185 644 444 739 5139(-16)
besselK_inc_ite(14.75, 19, 0, traceIBF = TRUE, epsilon = 0.95, nmax = 160)
## x = 17, y = 15, lambda = 1, exact value= 1.917 488 390 220 793 6555(-15)
besselK_inc_ite(17, 15, 1, traceIBF = TRUE, epsilon = 0.95, nmax = 160)
} 
