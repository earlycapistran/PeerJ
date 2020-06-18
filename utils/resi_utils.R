# =============================================================================
# Utilities for residual analysis 
# Supporting Information for Early-Capistrán et al. (2020), PeerJ
# earlycapistran@comunidad.unam.mx - April 2020
# =============================================================================

# Functions for residual analysis for gls and nls objects
# Developed by M.M. Early-Capistrán and H. Salinas.

# .............................................................................
# resi.analysis returns residual plots (normality, residuals vs. fitted, 
# autocorrelation), residual mean, Shapiro-Wilk normality test, Levene Test 
# for homogeneity of variance, and Run's test for randomness
# .............................................................................

# To run resi.analysis, you must have dplyr, car, and DescTools installed
resi.analysis  <- function(model) {
  resi <- residuals(model)
  resDf <- as.data.frame(resi) %>% select(resi)
  resDf$sign <- as.factor(ifelse(resDf < 0, "negative", "positive"))
  fit <- fitted(model)
  # Plot normality
  par(mfrow = c(2, 2))
  qqnorm(resDf$resi)
  qqline(resDf$resi)
  # Plot residual vs. fitted values
  res.plot <- plot(x = fit, y = resi,
       xlab = "Fitted values", ylab = "Residuals",
       main = "Residuals versus fitted values")
    # Plot residual autocorrelation
  plot(resDf$resi, c(resDf$resi[-1], +
                                  NA),
       xlab = "Residuals", ylab = "Lagged residuals",
       main = "Autocorrelation") 
  par(mfrow = c(1, 1)) 
  # Get residual mean
  mean <- mean(resDf$resi)
  # Run tests
  norm <- shapiro.test(resDf$resi)
  levene <- leveneTest(resDf$resi ~ sign,
             data = resDf)
  runs <- RunsTest(resDf$resi)
  result <- list(normality = norm, mean = mean, levene = levene, runs = runs)
  names(result) <- c("Residual Normality Test", "Residual Mean", "Levene's Test", "Runs Test")
  print(result)
}

#..............................................................................
# One-sample t-test (residual mean = 0) for gls models
# .............................................................................

# To run resi.t.test.gls you must have dplyr installed
resi.t.test.gls  <- function(model) {
  # Extract residuals
  resi <- residuals(model)
  resDf <- as.data.frame(resi) %>% select(resi)
  # Run t-test for mean = 0 ...................................................
  # Get residual standard deviation
  stDev <- sd(resDf$resi)
  mean <- mean(resDf$resi)
  # Get degrees of freedom
  degF <-summary(model)$dims$N 
  # Calculate t-value
  t.value  <- abs(mean/stDev)
  # Calculate p-value
  p.value <- dt(t.value, df=degF)
  # Print results of t-test 
  cat(paste("T-test", '\n', "t-value", t.value, ",", "p-value", p.value, '\n',
            "Null hypothesis: true mean is equal to 0",  '\n',
            "Alternative hypothesis: true mean is not equal to 0",'\n'))
}

#..............................................................................
# One-sample t-test (residual mean = 0) for nls models
# .............................................................................

# To run resi.t.test.nls you must have dplyr installed
resi.t.test.nls  <- function(model) {
  # Extract residuals
  resi <- residuals(model)
  resDf <- as.data.frame(resi) %>% select(resi)
  # Run t-test for mean = 0 ...................................................
  # Get residual standard deviation
  stDev <- sd(resDf$resi)
  mean <- mean(resDf$resi)
  # Get degrees of freedom
  degF <- summary(model)$df[2]
  # Calculate t-value
  t.value  <- abs(mean/stDev)
  # Calculate p-value
  p.value <- dt(t.value, df=degF)
  # Print results of t-test 
  t.result <- cat("T-test", '\n', "t-value", t.value, ",", "p-value", p.value, '\n',
                  "Null hypothesis: true mean is equal to 0",  '\n',
                  "Alternative hypothesis: true mean is not equal to 0",'\n'
  )
}