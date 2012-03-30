besselK_inc_clo <- function (x, z, lambda, lower = FALSE, 
  expon.scaled = FALSE) {
  is.wholenumber <- function(x, tol = .Machine$double.eps^0.5)
    abs(x - round(x)) < tol
  if (z >= 705 & expon.scaled == FALSE)
    warning("Underflow occurs, expon.scaled should be TRUE")
  if (z == 0)
    stop ("z should be positive")
  if ( x < 0)
    stop("x should be greater than zero")
  if (lambda < 0)
    stop("lambda should be positive")
  if (!is.wholenumber(lambda + 1/2))
    stop("lambda should be half of an odd positive integer")
  s <- lambda + 1/2
  Az <- sqrt(pi/(2 * z)) * exp(-z)
  if (expon.scaled == TRUE)
    Az <- sqrt(pi/(2 * z))
  r <- seq(0, s-1, length.out = s)
  evalbk <- rep(0, length(r))
  for(i in 1:length(evalbk)) {
    alpha <- lambda + r[i] + 1/2
    beta <- factorial(lambda - r[i]- 1/2)
    varphi <- seq(0, alpha - 1, length.out = alpha)
    binoCoef <- choose(alpha-1, varphi)
    binoVal <- numeric(length(binoCoef))
    for (j in 1:length(binoVal)) {
      binoVal[j] <- binoCoef[j] * x^(-varphi[j]) * 
      factorial(varphi[j]) 
    }
    if (lower== FALSE) {
      evalbk[i] <- 1/(beta*factorial(r[i]))*
        (1/(2 * z))^(r[i]) * x^(alpha - 1) * exp(-x) * sum(binoVal)
    } else {
      evalbk[i] <- 1/(beta*factorial(r[i]))*
        (1/(2 * z))^(r[i]) * (factorial(alpha - 1) - 
        x^(alpha - 1)*exp(-x) * sum(binoVal))
    }   
  }
  Az * sum(evalbk)
}

besselK_app_ser <- function(z, lambda, details = TRUE) { 
  lambda <- abs(lambda)
  maxNum <- 1.797693134862316e+308 
  minNum <- 2.220446049250313e-16 
  k <- 1  
  InitKterm <- 1  
  InitTerm <- 0 
  maxDenom <- 1
  conv <- 1
  maxTol <-  maxNum^(0.98)     
  largeVal <- function(z, lambda) { 
    if(z < 0) stop("z should be greater than zero")
    lambda <- abs(lambda)
    Xterm <- sqrt(pi/(2 * z)) * exp(-z)      
    beKex <- besselK(z, lambda)  
    Xterm <- sqrt(pi/(2 * z)) * exp(-z) 
    beKex <- besselK(z, lambda)  
    if (z >= 706) {
      Xterm <- sqrt(pi/(2 * z))
      beKex <- besselK(z, lambda, expon.scaled = TRUE) 
    }     
    while(maxDenom < maxTol) {  
      k <- k + 1 
      if (1/conv <= minNum^(0.98)) {
        zlarge <- 0
          break
      } else {
        zlarge <- 1
      }      
      Kterm <- 4 * lambda^2 - (2 * k - 3)^2 
      maxDenom <- factorial(k - 1) * (8 * z)^(k-1) 
      deno <- (factorial(k-1) * (8 * z)^(k-1))^(-1)   
      term1 <-  InitKterm * Kterm * deno  
      beK <- Xterm * sum(1, InitTerm, term1) 
      InitKterm <- Kterm * InitKterm 
      InitTerm <- term1 + InitTerm 
      conv <- abs((beKex - beK))  
    }
    return(zlarge)
  }
  is.wholenumber <- function(z, tol1 = minNum^(0.5))
  abs(z - round(z)) < tol1
  if ( z < 0)
    stop("z must be greater than zero")
  if ( z > maxNum^(0.98))
    stop("z should be less than maximum floating number")
  if (lambda > 90 & !is.wholenumber(lambda + 1/2))
    stop("Absolute value of lambda should be <= 90")
  if (z < 0)
    stop("z should be greater than zero")
  sigma <- 1
  Xterm <- sqrt(pi/(2 * z)) * exp(-z) 
  beKex <- besselK(z, lambda)  
  if (z < 705) {
    xVec <- seq(z, 720, by = 1)
    if(beKex <= minNum^(0.95)) sigma <- 1/beKex
    tol <- minNum^(0.95)  
    if(lambda >= 34 & lambda < 39) tol <- minNum^(0.7) 
    if(lambda >= 39 & lambda < 45) tol <- minNum^(0.5) 
    if(lambda >= 45 & lambda < 50) tol <- minNum^(0.3) 
    if(lambda >= 50 & lambda <= 90) tol <- minNum^(0.22) 
    beKex <- beKex * sigma
  }
  if (z >= 705) {
    if (!is.wholenumber(lambda + 1/2)) 
      cat("Approximated value is on expon.scale","\n") 
    Xterm <- sqrt(pi/(2 * z))
    beKex <- besselK(z, lambda, expon.scaled = TRUE) 
    tol <- minNum^(0.95) * sigma
    beKex <- beKex * sigma  
  } 
  k <- 1  
  InitKterm <- 1  
  InitTerm <- 0 
  conv <- 1
  if (is.wholenumber(lambda + 1/2)) {
    nterms <- lambda - 1/2 
    if (lambda == 1/2) nterms <- 1 
    beKex <- besselK(z, lambda, expon.scaled = FALSE) 
    if (z >= 705) {
          beKex <- besselK(z, lambda, expon.scaled = TRUE) 
          cat("Evaluation is on expon.scale", "\n")    
      }
    cat("lambda             :", lambda, "\n") 
    cat("Exact evaluation   :", beKex, "\n") 
    cat("Number of terms    :", lambda + 1/2, "\n")
  } else {
    if (z < 705) {
      zval <- largeVal(z, lambda)
      if (zval == 0) {
        for (i in 1:length(xVec)) {
          zLarge <- largeVal(xVec[i], lambda)
          if (identical(zLarge, 1)) {
            minVal <- xVec[i]
              break
          }                        
        }
          cat("For absolute value of lambda =", lambda, 
              ",z should be >=", minVal, sep = " ","\n")
      }             
    } else {
      zval <- 1
    } 
    if (zval == 1) {
      while (conv > tol) {  
      k <- k + 1 
      conv.old <- conv
      Kterm <- 4 * lambda^2 - (2 * k - 3)^2  
      deno <- (factorial(k-1) * (8*z)^(k-1))^(-1)  
      term1 <-  InitKterm * Kterm * deno  
      beK <- (Xterm * sum(1, InitTerm, term1)) * sigma
      InitKterm <- Kterm * InitKterm 
      InitTerm <- term1 + InitTerm 
      conv <- abs(beKex  - beK)
      }
      if (details == TRUE) {      
      cat("Approximated value         :", beK/sigma, "\n") 
      cat("Approximation error        :", conv, "\n")  
      cat("Number of terms required   :", k - 1, "\n")
      } else {
        approBK <- beK/sigma
        return(approBK)                     
      }
    }        
  }  
}
