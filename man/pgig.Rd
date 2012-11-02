\name{pgig}
\alias{pgig}

\title{Accurate Evaluation of Tail Probabilities of the Generalized
  Inverse Gaussian Distribution}

\description{
  Evaluates analytical formulae for distribution functions of the
    generalized inverse Gaussian distribution (GIG) by calling
    function \code{besselK_inc_err}.
}

\usage{
   pgig(q, lambda, chi, psi, lower.tail = TRUE, bit = 200)
}

\arguments{
  \item{q}{Quantile, \eqn{q > 0}.}
  \item{lambda}{Parameter, \eqn{\lambda = \pm(j + \frac{1}{2})}
    with \eqn{j = 0, 1, 2, \ldots}.}
  \item{chi}{Parameter chi > 0.}
  \item{psi}{Parameter psi > 0.}
  \item{lower.tail}{Logical. \eqn{P(W < w)} is
    returned if \code{TRUE}.}
  \item{bit}{Precision bit. A positive integer greater or equal 100.}	
}

\details{
 The GIG is given by
 \deqn{GIG(w|\lambda, \chi, \psi) =
 \frac{(\psi/\chi)^{\lambda/2}}{2 K_{\lambda}
 (\sqrt{\chi\psi})}\,e^{-\left(\chi w^{-1} \,+\,
 \psi w \right)/2}\, w^{\lambda - 1} \qquad w >0.}
 This distribution has been used in hydrology, reliability
 analysis, extreme events modelling in
 financial risk management, and as the mixing
 distribution to form the family of generalized hyperbolic
 distributions in statistics.
}

\references{

\cite{Barndorff-Nielsen, O. E (1977) Exponentially decreasing
  distributions for the logarithm of particle size.
  \emph{Proceedings of the Royal Society of London.
  Series {A}}, \bold{353}, 401--419.}

\cite{Olver, F.W.J., Lozier, D.W., Boisver, R.F. and Clark,
  C.W (2010) Handbook of Mathematical Functions.
  New York: National Institute of Standards and Technology,
  and Cambridge University Press.}

 \cite{Tran, T. T., Yee, W.T. and Tee, J.G (2012)
  Formulae for the Extended Laplace Integral and
  Their Statistical Applications. Working Paper.}

  \cite{Watson, G.N (1931) A Treatise on the Theory
  of Bessel Functions and Their Applications to Physics.
  London: MacMillan and Co.}

}

\note{
  This function allows for accurate evaluation of distribution
  functions (c.d.f and c.c.d.f) of the family of GIG distributions
  with \eqn{\lambda = \pm(j + \frac{1}{2})}. Currently,
  only c.d.f of inverse Gaussian
  distribution, \eqn{\lambda = -\frac{1}{2}}, is available.
}

\seealso{\code{CalIncLapInt}, \code{gamma_inc_err}}

\author{
  Thanh T. Tran \email{frmqa.package@gmail.com}
}

\examples{
## Accuracy tests
q <- 1
chi <- 3
psi <- 15
lambda <- 5/2
lowerTail <- sapply(lambda, function(w.)
  pgig(q, chi, psi, lambda = w., lower.tail = TRUE, 200))
upperTail <- sapply(lambda, function(w.)
  pgig(q, chi, psi, lambda = w., lower.tail = FALSE, 200))
## sum of two parts equals 1
(lowerTail + upperTail)

}

