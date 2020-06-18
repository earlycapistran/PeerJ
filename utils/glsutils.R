# =============================================================================
# Utilities for GLS objects
# Supporting Information for Early-Capistrán et al. (2020), PeerJ
# earlycapistran@comunidad.unam.mx - April 2020
# =============================================================================

# Functions for Durbin-Watson test and D^2 value for gls objects
# Developed by M.M. Early-Capistrán and H. Salinas.

#..............................................................................
# Durbin-Watson test for gls objects
# Modified from durbinWatsonTest in "car" package:
# https://rdrr.io/cran/car/src/R/durbinWatsonTest.R
# You must have the car package installed to run this function
#..............................................................................

durbinWatsonTest.gls <- function(model, data, max.lag=1, simulate=TRUE, reps=10000, 
                                 method=c("resample","normal"), 
                               alternative=c("two.sided", "positive", "negative"), ...){
  attach(data)
  method <- match.arg(method)
  alternative <- if (max.lag == 1) match.arg(alternative)
  else "two.sided"
  residuals <-as.vector(residuals(model))
  if (any(is.na(residuals))) stop ('residuals include missing values')
  n <- length(residuals)
  r <- dw <-rep(0, max.lag)
  den <- sum(residuals^2)
  for (lag in 1:max.lag){
    dw[lag] <- (sum((residuals[(lag+1):n] - residuals[1:(n-lag)])^2))/den
    r[lag] <- (sum(residuals[(lag+1):n]*residuals[1:(n-lag)]))/den
  }
  if (!simulate){
    result <- list(r=r, dw=dw)
    class(result) <- "durbinWatsonTest"
    result
  }
  else {
    S <- summary(model)$sigma
    X <- model.matrix(model)
    mu <- fitted(model)
    Y <- if (method == "resample") 
      matrix(sample(residuals, n*reps, replace=TRUE), n, reps) + matrix(mu, n, reps)
    else matrix(rnorm(n*reps, 0, S), n, reps) + matrix(mu, n, reps)
    E <- residuals(lm(Y ~ X - 1))
    DW <- apply(E, 2, car::durbinWatsonTest, max.lag=max.lag)
    if (max.lag == 1) DW <- rbind(DW)
    p <- rep(0, max.lag)
    if (alternative == 'two.sided'){
      for (lag in 1:max.lag) {
        p[lag] <- (sum(dw[lag] < DW[lag,]))/reps
        p[lag] <- 2*(min(p[lag], 1 - p[lag]))
      }
    }
    else if (alternative == 'positive'){
      for (lag in 1:max.lag) {
        p[lag] <- (sum(dw[lag] > DW[lag,]))/reps
      }
    }
    else {
      for (lag in 1:max.lag) {
        p[lag] <- (sum(dw[lag] < DW[lag,]))/reps
      }
    }
    detach(data)
    result <- list(r=r, dw=dw, p=p, alternative=alternative)
    class(result)<-"durbinWatsonTest"
    result
  }
}

print.durbinWatsonTest <- function(x, ...){
  max.lag <- length(x$dw)
  result <- if (is.null(x$p)) cbind(lag=1:max.lag,Autocorrelation=x$r, "D-W Statistic"=x$dw)
  else cbind(lag=1:max.lag,Autocorrelation = x$r, "D-W Statistic" = x$dw, 
             "p-value"= x$p)
  rownames(result) <- rep("", max.lag)
  print(result)
  cat(paste(" Alternative hypothesis: rho", if(max.lag > 1) "[lag]" else "",
            c(" != ", " > ", " < ")[which(x$alternative == c("two.sided", "positive", "negative"))],
            "0\n", sep=""))
  invisible(x)
}

#..............................................................................
# Calculate D^2 value for gls models
#..............................................................................

d.sq <- function(model) {
  null.deviance  <-  as.numeric(abs(((-2)*logLik(model))))
  res.deviance <- abs(sum(residuals(model)^2))
  d2 <- (null.deviance-res.deviance)/null.deviance
  return(d2)
}
