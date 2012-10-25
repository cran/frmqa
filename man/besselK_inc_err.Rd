\name{besselK_inc_err}
\alias{besselK_inc_err}

\title{Calculation of the Incomplete BesselK Functions in Terms of the Complementary Error Functions}

\description{
  Calculates incomplete BesselK functions by evaluating explicit 
  expressions for the lower and upper incomplete BesselK in terms 
  of the complementary error function by calling \code{CalIncLapInt}.
}

\usage{
   besselK_inc_err(x, z, lambda, bit, lower = FALSE)
}

\arguments{
  \item{x}{Argument, \code{x > 0}.}
  \item{z}{Argument, \code{z > 0}.}
  \item{lambda}{Argument, \eqn{\lambda = \pm(j + \frac{1}{2})} 
    with \eqn{j = 0, 1, 2, \ldots}.}
  \item{lower}{Logical. Lower incomplete Bessel function 
    is calculated if \code{TRUE}.}
  \item{bit}{Precision bit. A positive integer greater or equal 100.}	
}

\details{
  One of the integral representations of the lower incomplete 
  BesselK is given by 
  \deqn{\widehat K_{\lambda}(z, x) = 
  \frac{1}{\left(2z\right)^{\lambda}}\,\int_0^x\, 
  e^{-\left\{z^2\,\xi^2 \,+\, 1/(4\, \xi^2)\right\}}\,
  \xi^{-2\lambda -1}\, d\xi,} which appears in the distribution 
  function of the generalized inverse Gaussian distribution, 
  see Barndorff-Nielsen(1977).       
}

\references{

\cite{Barndorff-Nielsen, O. E (1977) Exponentially decreasing 
  distributions for the logarithm of particle size. \emph{Proceedings 
  of the Royal Society of London. Series {A}}, \bold{353}, 401--419.}

\cite{Olver, F.W.J., Lozier, D.W., Boisver, R.F. and Clark,
  C.W (2010) Handbook of Mathematical Functions.
  New York: National Institute of Standards and Technology,
  and Cambridge University Press.}

 \cite{Tran, T. T., Yee, W.T. and Tee, J.G (2012) Formulae for the Extended
      Laplace Integral and Their Statistical Applications. Working Paper.}

  \cite{Watson, G.N (1931) A Treatise on the Theory of Bessel 
  Functions and Their Applications to Physics. London: MacMillan and Co.}

}

\note{
  Currently, analytical formulae for the incomplete 
  BesselK functions are not available for any value 
  of \code{lambda}. 
}

\seealso{
  \code{besselK_inc_clo}, \code{gamma_inc_clo}, \code{pgig} 
}

\author{
  Thanh T. Tran \email{frmqa.package@gmail.com}
}

\examples{
  ## Accuracy tests
  x <- 2
  z <- 5
  lambda <- -c(1/2, 5/2)
  lower <- sapply(lambda, function(w.)
    besselK_inc_err(x, z, lambda = w., 200, lower = TRUE))
  upper <- sapply(lambda, function(w.)
    besselK_inc_err(x, z, lambda = w., 200, lower = FALSE)) 
  ## sum of two parts   
  (lower + upper)
  ## equals the whole function
  (besselK(z, nu = lambda))
  
}
