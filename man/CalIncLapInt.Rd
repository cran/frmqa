\name{CalIncLapInt}
\alias{CalIncLapInt}

\title{Evaluation of Analytical Formulae for the 
  Incomplete Laplace Integral}

\description{
    Evaluates analytical formulae for the 
    lower and upper incomplete Laplace integral in terms of the 
    complementary error function. 
}
\usage{
   CalIncLapInt(lambda, a = 1, b = 1, x = 1, lower = TRUE, bit = 200)   
}

\arguments{
    \item{lambda}{Order \eqn{\lambda = \pm(j + \frac{1}{2})} 
        with \eqn{j = 0, 1, 2, \ldots}.}
    \item{a}{Argument a > 0.}
    \item{b}{Argument b >= 0.}
    \item{x}{Limit of the incomplete Laplace integral.}
    \item{lower}{Logical. Lower incomplete Laplace integral 
      is returned if \code{TRUE}.}
    \item{bit}{Precision bit. A positive integer greater 
    or equal 100.}	
}

\details{The lower and upper extended Laplace integrals 
  are given by \deqn{\widehat L_{\lambda}(x, a, b) =  \int_0^{x} 
  e^{-(a \xi^2 + b/\xi^2)}\xi^{-2\lambda - 1}\; d\xi ,}  and 
  \deqn{\widetilde L_{\lambda}(x, a, b) =  \int_x^{\infty} 
  e^{-(a \xi^2 + b/\xi^2)}\xi^{-2\lambda - 1}\; d\xi} respectively. 
  Calculation is performed using multiple precision 
  floating-point reliably or MPFR-numbers instead of the 
  default floating-point number in \bold{R}, which ensure 
  accuracy is at least 100 bit.
}

\references{
  \cite{Hankin, R.K.S (2006) {A}dditive integer partitions in {R}. 
  Journal of Statistical Software, Code Snippets, \bold{16}.}

\cite{Maechler, M Rmpfr: R MPFR - Multiple Precision Floating-Point 
  Reliable. R package version 0.5-0, \url{http://CRAN.R-project.
  org/package=Rmpfr}.}

\cite{Olver, F.W.J., Lozier, D.W., Boisver, R.F. and Clark,
  C.W (2010) Handbook of Mathematical Functions. 
  New York: National Institute of Standards and Technology, 
  and Cambridge University Press.}
  
\cite{Tran, T. T., Yee, W.T. and Tee, J.G (2012) Formulae for 
    the Extended Laplace Integral and Their Statistical 
    Applications. Working Paper.}
      
  \cite{Watson, G.N (1931) A Treatise on the Theory of 
    Bessel Functions and Their Applications to 
    Physics. London: MacMillan and Co.}
}

\author{
  Thanh T. Tran \email{frmqa.package@gmail.com}
}        

\seealso{
   \code{besselK_inc_err}, \code{gamma_inc_err}, \code{pgig} 
}
