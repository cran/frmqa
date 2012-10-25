\name{BesselK_inc_erfc}
\alias{BesselK_inc_erfc}

\title{Symbolic Calculation of Lower and Upper Incomplete Bessel Function}

\description{
    This function evaluates lower and upper incomplete Bessel function, see details for the integral representation in this case, by calling function \code{CalIncLapInt}.
}

\usage{
    BesselK_inc_erfc(lambda, z, x, lower = TRUE, bit = 200)
}

\arguments{
    \item{lambda}{Order. Equals \eqn{\pm(j + \frac{1}{2})} with \eqn{j = 0, 1, 2, \ldots}.}
    \item{z}{Argument z > 0.}
    \item{x}{Argument x > 0. Limit of the incomplete Bessel function.}
    \item{lower}{Logical. Lower incomplete Bessel function is returned if TRUE.}
    \item{bit}{Precision bit. A positive integer greater or equal 100.}	
}

\details{The lower and upper incomplete Bessel functions are given by
  \deqn{\widehat K_{\lambda}(z, x) = \frac{1}{(2\,z)^{\lambda}} \int_0^{x} e^{-(a \xi^2 + b/\xi^2)}\xi^{-2\lambda - 1}\;
	 d\xi ,} and \deqn{\widetilde K_{\lambda}(z, x) = \frac{1}{(2\,z)^{\lambda}} \int_x^{\infty} e^{-(a \xi^2 + b/\xi^2)}\xi^{-2\lambda - 1}\;	 d\xi} respectively. These incomplete Bessel functions appear in distribution functions of the family of the Generalized Inverse Gaussian (GIG) and Generalized Hyperbolic distributions.}

\references{

\cite{Olver, F.W.J., Lozier, D.W., Boisver, R.F. and Clark,
  C.W. (2010) Handbook of Mathematical Functions.
  New York: National Institute of Standards and Technology,
  and Cambridge University Press.}

  \cite{Tran, T. T. (2011) Some Problems Concerning the Generalized Hyperbolic and Related Distributions. Ph.D Thesis.
        University of Auckland, New Zealand.}

  \cite{Tran, T. T., Yee, W.T. and Tee, J.G.(2012) Formulae for the Extended
      Laplace Integral and Their Statistical Applications. Working Paper.}

  \cite{Watson, G.N (1931) A Treatise on the Theory of Bessel Functions and Their Applications to Physics.
  London: MacMillan and Co.}
}

\note{
    Analytical formulae for evaluating \eqn{\widehat K_{\lambda}(x, a, b)} and \eqn{\widetilde K_{\lambda}(x, a, b)} have not been documented in the literature for any value of \code{lambda}. Even Maple 15 cannot symbolically evaluate them in the simplest cases where \eqn{\lambda = \pm\frac{1}{2}}. Evaluation of \eqn{\widehat K_{\lambda}(x, a, b)} and \eqn{\widetilde K_{\lambda}(x, a, b)} using \code{BesselK_inc_erfc} allows for exact calculation of tail probabilities of so-called \bold{exceedances} in modelling
    extreme events using the GIG. Note that \eqn{\widehat K_{\lambda}(z, x)} and \eqn{\widetilde K_{\lambda}(z, x)} converge to \eqn{K_{\lambda}(z)} when \eqn{x = \infty} or \eqn{x = 0} respectively.
}

\seealso{\code{CalIncLapInt}, \code{besselK_inc_clo}}

\author{
  Thanh T. Tran \email{frmqa.package@gmail.com}
}


