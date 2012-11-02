\name{gamma_inc_err}
\alias{gamma_inc_err}

\title{Accurate Calculation of the Incomplete Gamma Functions 
  Using Analytical Formulae}

\description{
  Evaluates explicit formulae for the lower and upper 
  incomplete gamma functions in terms of complementary 
  error function by calling \code{CalIncLapInt}.
}

\usage{
   gamma_inc_err(x, lambda, bit, lower = FALSE)
}

\arguments{
  \item{x}{Argument, \code{x > 0}.}
  \item{lambda}{Argument, \eqn{\lambda = \pm(j + \frac{1}{2})} 
  with \eqn{j = 0, 1, 2, \ldots}.} 
  \item{lower}{Logical. Lower incomplete gamma function is 
    calculated if \code{TRUE}.}
  \item{bit}{Precision bit. A positive integer greater or equal 100.}	
}

\details{
  The lower incomplete gamma function is given by
  \deqn{\gamma(x, \lambda) = \int_0^x e^{-t}\,t^{\lambda - 1}\, dt.}
}

\references{

\cite{Olver, F.W.J., Lozier, D.W., Boisver, R.F. and Clark,
  C.W (2010) Handbook of Mathematical Functions.
  New York: National Institute of Standards and Technology,
  and Cambridge University Press.}

  \cite{Tran, T. T (2011) Some Problems Concerning the 
    Generalized Hyperbolic and Related Distributions. Ph.D Thesis.
    The University of Auckland, New Zealand.}

  \cite{Tran, T. T., Yee, W.T. and Tee, J.G (2012) Formulae 
    for the Extended Laplace Integral and Their Statistical 
    Applications. Working Paper.}

  \cite{Watson, G.N (1931) A Treatise on the Theory of Bessel 
  Functions and Their Applications to Physics.
  London: MacMillan and Co.}

}

\note{
  This function evaluates formulae in terms of complementary 
  error function for \eqn{\gamma(x, \lambda)} and its 
  upper counterpart when \eqn{\lambda = \pm(j + \frac{1}{2})}. 
  Currently, such formulae are only available when 
  \eqn{\lambda = \pm\frac{1}{2}}.
}

\seealso{
  \code{CalIncLapInt}, \code{besselK_inc_clo}, \code{pgig}
}

\author{
  Thanh T. Tran \email{frmqa.package@gmail.com}
}

\examples{
  ## Accuracy tests
  x <- 3
  lambda <- 3/2
  lower <- sapply(lambda, function(w.)
    gamma_inc_err(x, lambda = w., 200, lower = TRUE))
  upper <- sapply(lambda, function(w.)
    gamma_inc_err(x, lambda = w., 200, lower = FALSE))  
 ## sum of two parts   
  (lower + upper)
 ## equals the whole function 
  (gamma(lambda))
  
}
