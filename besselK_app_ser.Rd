\name{besselK_app_ser}
\alias{besselK_app_ser}

\title{Asymptotic Approximations of the Modified Bessel K Function}

\description{
    Approximates modified Bessel function of the third 
    kind \eqn{K_{\lambda}(z)} using asymptotic series.
}

\usage{
    besselK_app_ser(z, lambda, details = TRUE)
}

\arguments{
    \item{z}{Numeric, \code{z > 0}. }
    \item{lambda}{Real valued in the interval \eqn{[-90,90]}.}
    \item{details}{Logical. If \code{TRUE} then approximation 
    details are returned.}
}

\details{ 
  The asymptotic series of \eqn{K_{\lambda}(z)} can 
  be written as \deqn{K_{\lambda}(z) \sim  \sqrt{\frac{\pi}{2 z}}\, 
  e^{-z} \Bigg\{1 + \frac{4\lambda^2 - 1^2}{1!\, 8 z} + 
  \frac{\left(4 \lambda^2 - 1^2\right)\left(4 \lambda^2 - 
  3^2 \right)}{2!\,(8 z)^2} + \cdots \Bigg\}.} For a given 
  value of \eqn{\lambda}, this function checks if the 
  entered value of \code{z} is sufficiently large. Otherwise it 
  approximates this value. The function returns: an accurate 
  approximated value of \eqn{K_{\lambda}(z)} (up to 19 dp) and a 
  finite number of terms in the summation. 
}

\references{
  \cite{Abramowitz, M., Stegun, I. A., (1970) 
  Handbook of Mathematical Functions with Formulas, 
  Graphs, and Mathematical Tables. New York: Dover.}
    
  \cite{Olver, F.W.J., Lozier, D.W., Boisver, R.F. and Clark,
  C.W. (2010) Handbook of Mathematical Functions. 
  New York: National Institute of Standards and Technology, 
  and Cambridge University Press.}
  
  \cite{Tran, T. T. (2011) Some Problems Concerning the Generalized Hyperbolic 
  and Related Distributions. Ph.D Thesis. University of Auckland, New Zealand.}
}  

\author{
  Thanh T. Tran \email{frmqa.package@gmail.com}
}   

\note{
  By using \code{besselK_app_ser}  
  formulae for mean, variance, skewness 
  and kurtosis of GH and GIG can be expressed in closed-form if the shape parameter
  is sufficiently larege. If \code{z >=705} then approximated value is on an 
  exponential scale to avoid underflow.}

\examples{
  besselK_app_ser(73, 9.12, details = TRUE)
  besselK_app_ser(1, 5, details = TRUE)
}
