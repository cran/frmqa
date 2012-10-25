\name{besselK_inc_clo}
\alias{besselK_inc_clo}
\title{Exact Calculation of the Incomplete BesselK Function}

\description{
  Calculates upper and lower incomplete functions of the 
  modified Bessel function of the third kind
  \eqn{K_{\lambda}(z)}, see details, using closed-form formulae.
}

\usage{
  besselK_inc_clo(x, z, lambda, lower = FALSE, expon.scaled = FALSE)
}

\arguments{
  \item{x}{Limit of the integration, \code{x > 0}.}
  \item{z}{Argument of the function, \code{z > 0}.}
  \item{lambda}{Order, \eqn{\lambda = (j + \frac{1}{2})} 
    with \eqn{j = 0, 1, 2, ...}.}
  \item{lower}{Logical. If \code{TRUE} then the lower incomplete 
    BesselK is calculated.}
  \item{expon.scaled}{Logical. If \code{TRUE} then the result 
    is on an exponential scale.}
}

\details{
  One of the integral representations of \eqn{K_{\lambda}(z)} 
  is given by \deqn{K_{\lambda}(z) = \sqrt{\frac{\pi}{2z}}
  \frac{1}{\Gamma(\lambda + \frac{1}{2})}\,e^{-z}\int_0^{\infty}
  e^{-\xi}\xi^{\lambda - 1/2}\left(1+\frac{\xi}{2 z}
  \right)^{\lambda -1/2}\,d\xi,} \code{besselK_inc_clo} evaluates 
  closed-form formulae, which we derived to compute this 
  integral, in the \eqn{(0, x)} and \eqn{(x, \infty)} intervals
  for the so-called lower and upper incomplete Bessel 
  function respectively. ``Exact" evaluation of the integral 
  in these intervals can also be obtained by 
  numerical integration using software such 
  as Maple \url{www.maple.com}.
}
\references{
   \cite{Olver, F.W.J., Lozier, D.W., Boisver, R.F. and Clark,
    C.W (2010) Handbook of Mathematical Functions. 
    New York: National Institute of Standards and Technology, 
    and Cambridge University Press.}
  
  \cite{Watson, G.N (1931) A Treatise on the Theory of 
    Bessel Functions and Their Applications to Physics. 
    London: MacMillan and Co.}
}

\author{Thanh T. Tran \email{frmqa.package@gmail.com}}

\seealso{
  \code{besselK_app_ser}, \code{besselK_inc_erfc}
}

\examples{
  options(digits = 15)	
  ## For x = 5, z = 8, lambda = 15/2 Maple 15 gives exact value of the
  ## lower incomplete Bessel function 0.997 761 151 460 5189(-4)
  besselK_inc_clo(5, 8, 15/2, lower = TRUE, expon.scaled = FALSE)
  ## For x = 21, z = 8, lambda = 21/2 Maple 15 give exact value of the
  ## upper incomplete Bessel function 0.704 812 324 921 884 3938(-2)
  besselK_inc_clo(21, 8, 21/2, lower = FALSE, expon.scaled = FALSE)
}
