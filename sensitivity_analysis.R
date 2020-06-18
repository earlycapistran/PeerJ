# =============================================================================
# Sensitivity analysis
# Supporting Information for Early-Capistrán et al. (2020), PeerJ
# earlycapistran@comunidad.unam.mx - April 2020
# =============================================================================

# Based on Soetaert, K., & Petzoldt, T. (2010). Inverse Modelling, Sensitivity 
# and Monte Carlo Analysis in R Using Package FME. Journal of Statistical 
# Software, 33(3). https://doi.org/10.18637/jss.v033.i03

# =============================================================================
#  Install and load libraries
# =============================================================================
# Check if required libraries are installed and install if necessary ..........
packages <- c("FME", "ggplot2", "ggthemes")
if (length(setdiff(packages, rownames(installed.packages()))) > 0) {
  install.packages(setdiff(packages, rownames(installed.packages())))  
}

# Load libraries ..............................................................
library(FME)
library(ggplot2)
library(ggthemes)

# .............................................................................
# NOTE:
# If you did not open this script directly from ".../quantifying_lek_data_code",
# please define the working directory to ".../quantifying_lek_data_code":
#
# setwd(".../quantifying_lek_data_code")
#
# .............................................................................

# =============================================================================
#  Load and prepare data
# =============================================================================

# Load dataset from csv file ..................................................
st_data = read.csv("data/standardised_dataset.csv", header = TRUE)

# Remove rows with missing values .............................................
st_data = st_data[complete.cases(st_data),]

# Input observed values as data frame .........................................
obs <- data.frame(x=st_data$yearSerial,  # year (serialised)
                  y=st_data$stCpue)      # green turtle CPUE

# =============================================================================
#  Fitting the model to the data
# =============================================================================

# Define exponential model ....................................................
# Returns a data.frame, with elements x and y
nlrModel <- function(p, x) return(
  data.frame(x = x, y=p[1]*exp(p[2]*x)
  ))


# Estimate the deviances of the model vs. the data with the function
# "residuals"..................................................................
residuals  <- function(p) (obs$y - nlrModel(p, obs$x)$y)

# Inuut "residuals" to modFit, which fits the model to the observed data .....
# Starting values from NLR with best fit in Early Capistrán et al. (2020)
fit.model <- modFit(f = residuals, p = c(0.189111E+02 ,-.264181E+00)) 

# Estimate and print summary of fit ...........................................
sFit <- summary(fit.model)
sFit

# Plot residual sum of squares, residual, and best fit ........................
x      <-0:35
par(mfrow = c(2, 2))
plot(fit.model, mfrow = NULL)
plot(obs, pch = 16, cex = 2, xlim = c(0, 35), ylim = c(0, 20),
     xlab = "Year", ylab = "CPUE (turtles/night)", main = "Best fit")
lines(nlrModel(fit.model$par, x))
par(mfrow = c(1, 1))

# =============================================================================
#  MCMC analysis
# =============================================================================

# For the initial model variance (var0) we use the residual mean squares 
# returned by the summary function, with equal weight to prior and modeled mean 
# squares (wvar0=1)

CoVar <- sFit$cov.scaled 
s2prior <- sFit$modVariance

# Run MCMC at 5000 iterations .................................................
MCMC <- modMCMC(f = residuals, p = fit.model$par, jump = CoVar, niter = 5000,
                var0 = s2prior, wvar0 = 1)

# Update covariance to increase acceptance ....................................
MCMC <- modMCMC(f = residuals, p = fit.model$par, 
                jump = CoVar, niter = 5000,
                ntrydr = 3, var0 = s2prior, 
                wvar0 = 1, updatecov = 100)

MCMC$count

#Plot convergence .............................................................
plot(MCMC, Full = TRUE)

# Plot posterior distribution of parameters, sum of squares and model error 
# standard deviation ..........................................................
hist(MCMC, Full = TRUE, col = "darkblue")

# Plot pairs to visualise the relationship between the two parameters .........
pairs(MCMC)

# Caltulate  parameter  correlation  and  covariances  from  the  MCMC........
cor(MCMC$pars)
cov(MCMC$pars)
# These results can be compared to results obtained by the fitting algorithm:
sFit
sFit$cov.scaled

# Estimate posterior predictive distribution of the model .....................
# Note: this estimate only assumes parameter uncertainty
sR<-sensRange(parInput=MCMC$pars,func=nlrModel,x=1:35)

# Plot posterior predictive distribution ......................................
plot(summary(sR), quant = TRUE)
# Add observed data points
points(obs)

# Estimate measurement error by adding randomly distributed noise to the 
# model predictions produced by the parameters from the MCMC chain.............
# Extract parameter sets used to produce output in "sR"
pSet <- attributes(sR)$pset
# Add randomly distributed noise
nOut  <- nrow(sR)
sR2   <- sR
iVar  <- 3:ncol(sR)
error <- rnorm(nOut, mean = 0, sd = sqrt(MCMC$sig[pSet]))
sR2[,iVar] <- sR2[ ,iVar] + error

# Plot posterior predictive distribution with randomly distributed noise.......
plot(summary(sR2),quant=TRUE)
# Add observed data points
points(obs)

# Make custom plot with ggplot2 ...............................................
g2 <- ggplot(summary(sR2)) +
  geom_line(aes(x, Mean), linetype = 1) +
  geom_ribbon(aes(x, ymin = q05, ymax = q25), fill = "skyblue4", alpha = .5) +
  geom_ribbon(aes(x, ymin = q75, ymax = q95), fill = "skyblue4", alpha = .5) +
  geom_ribbon(aes(x, ymin = q25, ymax = q50), fill = "skyblue1", alpha = .5) +
  geom_ribbon(aes(x, ymin = q50, ymax = q75), fill = "skyblue1", alpha = .5) +
  geom_point(data=obs, aes(x=x, y=y)) +
  labs(y = "Mean CPUE (turtles/night)", x = "Year") +
  theme_hc() + 
  scale_color_economist()
g2
