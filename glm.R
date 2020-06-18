# =============================================================================
# GLMs with residual correlation structure
# Supporting Information for Early-Capistrán et al. (2020), PeerJ
# earlycapistran@comunidad.unam.mx - April 2020
# =============================================================================

# . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . 
# Based on:
# Zuur, A. F. (Ed.). (2009). Mixed effects models and extensions in 
# ecology with R. Springer.
# . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . 

# =============================================================================
#  Install and load libraries and scripts
# =============================================================================

# Check if required libraries are installed and install if necessary...........
packages <- c("nlme", "dplyr", "car", "DescTools")
if (length(setdiff(packages, rownames(installed.packages()))) > 0) {
  install.packages(setdiff(packages, rownames(installed.packages())))  
  }

library(nlme)
library(dplyr)
library(car)
library(DescTools)

# .............................................................................
# NOTE:
# If you did not open this script directly from ".../quantifying_lek_data_code",
# please define the working directory to ".../quantifying_lek_data_code":
#
# setwd(".../quantifying_lek_data_code")
# .............................................................................

source("utils/resi_utils.R")
source("utils/gls_utils.R")

# =============================================================================
#  Load and prepare data
# =============================================================================

# Load raw database from csv file .............................................
raw_data = read.csv("data/raw_dataset.csv", header = TRUE)

# Convert fisherCode to characters to provide unique values for each data point
raw_data$fisherCode <- as.character(raw_data$fisherCode)

# Check normality of response variable (raw CPUE) ............................
shapiro.test(raw_data$rawCpue)
# Check normality of log-transformed raw CPUE ................................
shapiro.test(raw_data$logRawCpue)

# =============================================================================
# Model 1
# =============================================================================

# Set up a model with auto-regressive correlation and temporal order of the 
# data specified by "year (serialised)". The variable "fisherCode" is 
# included as a grouping factor to assure unique values 
# within groups required for corAR1 objects.

cpue.mod.1 <- gls(logRawCpue ~ yearSerial + experience + vesselCapacity
                        -1,
                        na.action = na.omit, data = raw_data, 
                        correlation = corAR1(form =~ yearSerial |fisherCode))

summary(cpue.mod.1)
# Calculate D^2 ................................................................
d.sq(cpue.mod.1)

# -------------------------------------------------------------------------------
# Residual analysis
# -------------------------------------------------------------------------------

# Returns normality plot, plot of residual vs. fitted values, residual mean, 
# Shapiro-Wilk normality test, Levene Test for homogeneity of variance, and
# Run's test for randomness ....................................................
resi.analysis(cpue.mod.1)

# Test residual mean = 0 .......................................................
resi.t.test.gls(cpue.mod.1)

# Durbin-Watson test for residual independence ................................
durbinWatsonTest.gls(model = cpue.mod.1, data=raw_data)

# =============================================================================
# Model 2
# =============================================================================

cpue.mod.2 <- gls(logRawCpue ~ yearSerial + gear + totalNetLength + numNets
                + experience -1,
          na.action = na.omit, data = raw_data, 
          correlation = corAR1(form =~ yearSerial|fisherCode))

summary(cpue.mod.2)
d.sq(cpue.mod.2)

# -------------------------------------------------------------------------------
# Residual analysis
# -------------------------------------------------------------------------------

resi.analysis(cpue.mod.2)
resi.t.test.gls(cpue.mod.2)
durbinWatsonTest.gls(model = cpue.mod.2, data = raw_data)

# =============================================================================
# Model 3
# =============================================================================

cpue.mod.3 <- gls(logRawCpue ~ yearSerial + gear + netLength -1,
                 na.action = na.omit, data = raw_data, 
                 correlation = corAR1(form =~ yearSerial|fisherCode))

summary(cpue.mod.3)
d.sq(cpue.mod.3)

# -------------------------------------------------------------------------------
# Residual analysis
# -------------------------------------------------------------------------------
resi.analysis(cpue.mod.3)
resi.t.test.gls(cpue.mod.3)
durbinWatsonTest.gls(model = cpue.mod.3, data=raw_data)
