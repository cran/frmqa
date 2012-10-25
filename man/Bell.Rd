\name{Bell}
\alias{Bell}
\title{Bell Coefficients Calculation}

\description{
  This function calculates Bell coefficients as described 
  on p.134--139 of Comtet (1974) for calculating high order 
  derivative of composite function \eqn{f(g(x))} wrt \eqn{x}, say. 	
}

\details{
  \code{Bell} calls the functions in package \bold{partititions}, 
  see \url{http://cran.r-project.org/web/packages/partitions/partitions.pdf}. 
}

\arguments{
  \item{n}{Derivative order \code{n > 1}. }
}

\usage{
  Bell(n)
}

\note{
   Combinatorial computation intensifies as \code{n} increases.
}

\references{
  \cite{Comtet, L (1974) Advanced combinatorics: {T}he art of 
    finite and infinite expansions. Boston: D. Reidel Pub.}  
  
  \cite{Hankin, R.K.S (2006) {A}dditive integer partitions 
    in {R}. \emph{Journal of Statistical Software, Code Snippets}, \bold{16}.}
   
  \cite{Tran, T. T., Yee, W.T. and Tee, J.G (2012) Formulae for the 
    Extended Laplace Integral and Their Statistical 
    Applications. Working Paper.}    
}

\author{
  Thanh T. Tran \email{frmqa.package@gmail.com}
}
  \seealso{\code{derivErfc}
}
\examples{
  Bell(5)
}
