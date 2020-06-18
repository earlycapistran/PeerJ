# =============================================================================
# Nonlinear Regression
# Supporting Information for Early-Capistrán et al. (2020), PeerJ
# earlycapistran@comunidad.unam.mx - April 2020
# =============================================================================

# =============================================================================
#  Install and load libraries
# =============================================================================

# Check if required libraries are installed and install if necessary ..........
packages <- c("easynls", "nlstools", "dplyr", "ggplot2", "ggthemes", "car",
              "DescTools")
if (length(setdiff(packages, rownames(installed.packages()))) > 0) {
  install.packages(setdiff(packages, rownames(installed.packages())))  
}

# Load libraries ..............................................................
library(easynls)
library(nlstools)
library(car)
library(DescTools)
library(ggplot2)
library(ggthemes)
library(dplyr)

# .............................................................................
# NOTE:
# If you did not open this script directly from ".../quantifying_lek_data_code",
# please define the working directory to ".../quantifying_lek_data_code":
#
# setwd(".../quantifying_lek_data_code/")
# .............................................................................

source("utils/resi_utils.R")

# =============================================================================
#  Load and prepare data
# =============================================================================

# Load dataset from csv file ..................................................
st_data = read.csv("data/standardised_dataset.csv", header = TRUE)
# Standardised CPUE data can be found in ~/datasets/standardised_dataset.csv

# Select columns with X (serialYear) and Y (cpue) variables ...................
st_data <- select(st_data, yearSerial, stCpue)


# Remove rows with missing values .............................................
st_data = st_data[complete.cases(st_data),]

# =============================================================================
# Fitting the exponential model "y~a*exp(b*x)" to standardised CPUE values
# Starting values from LabFit 7.2.49 (http://labfit.net/)
# =============================================================================

# Define the model with nlstools ..............................................
exp.model <- as.formula(stCpue~a*exp(b*yearSerial))

# Define starting values and run model ........................................
cpue.exp <- nls(exp.model, start = list(a=0.189111E+02, b=-.264181E+00), 
                data = st_data, model=TRUE)

# Plot nonlinear regression ...................................................
plotfit(cpue.exp, smooth = TRUE, main="Exponential decay model",
        xlab="Year (serialised)", ylab="Mean standardised CPUE (turtles/night)")

# View results ................................................................
overview(cpue.exp)

# easynls provides an R^2 value ...............................................
nlsfit(st_data, model=6, start=c(0.189111E+02 ,-.264181E+00))

# =============================================================================
# Residual analysis                       
# =============================================================================

# Residual analysis function ..................................................
# Returns normality plot, plot of residual vs. fitted values, autocorrelation
# plot, residual mean, Shapiro-Wilk normality test, Levene Test for 
# homogeneity of variance, and Run's test for randomness 
resi.analysis(cpue.exp)

# Test residual mean = 0 ......................................................
resi.t.test.nls(cpue.exp)

# Evaluate residual autocorrelation with autocorrelation plot and correlation
# test (residuals vs. lagged residuals) .......................................
# Create a data frame with residuals and lagged residuals
resi <- residuals(cpue.exp) 
resDf <- as.data.frame(resi) %>% select(resi)
resDf$resi_lag <- c(resDf$resi[-1], NA) # Append NA to vector of lagged  
                                        # residuals to match vector lengths
# Generate autocorrelation plot of residuals with linear trend line
ggplot(aes(x=resi, y=resi_lag), data=resDf) + geom_point(size=2, na.rm=TRUE) + 
  geom_smooth(method = "lm", se=FALSE, na.rm=TRUE) +
  labs(y = "Lagged residuals", x = "Residuals")

# Run a linear regression of residuals vs. lagged residuals
lag_lm <- lm(resDf$resi ~ resDf$resi_lag)
summary(lag_lm)

# Run Pearson correlation test for residuals vs. lagged residuals
cor.test(resDf$resi, resDf$resi_lag,  method = "pearson")

# =============================================================================
# Custom NLR plot with ggplot                  
# =============================================================================

# Plot with ggplot ............................................................
custom_plot <- ggplot(st_data, aes(x = yearSerial, y = stCpue), ylim=20) +
  geom_point(size=2.5, colour = "steelblue", na.rm=TRUE) +
  geom_smooth(method="nls", formula = y ~ a*exp(b*x), na.rm=TRUE,
              method.args = list(start = c(a=0.189111E+02, b=-.264181E+00)), 
              se = F, colour = "grey35") + 
  labs(y = "Mean standardised CPUE (turtles/night)", x = "Year") +
  theme_hc()

custom_plot
