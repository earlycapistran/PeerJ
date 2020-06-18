# =============================================================================
# Lin Concordance Correlation Coefficient (Lin CCC)
# Supporting Information for Early-Capistrán et al. (2020), PeerJ
# earlycapistran@comunidad.unam.mx - February 2020
# =============================================================================

# =============================================================================
# Load libraries
# =============================================================================

# Check if required libraries are installed and install
# if necessary
packages <- c("DescTools", "ggplot2", "ggthemes", "dplyr")
if (length(setdiff(packages, rownames(installed.packages()))) > 0) {
  install.packages(setdiff(packages, rownames(installed.packages())))  
}

library(DescTools)
library(ggplot2)
library(ggthemes)
library(dplyr)

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

# Load standardised CPUE and fisheries landing data ...........................
data = read.csv("data/comp_dataset.csv", header=TRUE)

# Remove rows with missing values .............................................
comp_data = data[complete.cases(data),]
#attach(comp_data)

# Standardise values to z-scores ..............................................
# Define variables
cpue <- comp_data$stCpue
landing <- comp_data$totalLanding
# Define z-scores and add to data frame
comp_data[,'zStCpue']<-(cpue-mean(cpue))/sd(cpue)
comp_data[,'zLanding']<-(landing-mean(landing))/sd(landing)

# =============================================================================
#  Plot z-scores for exploratory visual evaluation
# =============================================================================

# Plot z-scores for standardised CPUE
p1 <- ggplot(comp_data, aes(x = yearSerial, y = zStCpue)) +
  geom_point(size=2.5, colour = "steelblue") +
  labs(x="Year", y="Standardised CPUE (z-score)") +
  theme_hc()+ scale_colour_hc()

#Plot z-scores of CPUE and total annual landings
p2 <- p1 +  geom_point(aes(x = yearSerial, y = zLanding), 
                       size=2.5, color="tomato", shape=1) +
  labs(x="Year (serialised)", y="z-scores for CPUE and annual landings") 
p2


# =============================================================================
#  Run Lin CCC
# =============================================================================

# Define variables ............................................................
zStCpue = comp_data$zStCpue
zLanding = comp_data$zLanding

# Run Lin CCC .................................................................
lin.ccc <- CCC(zStCpue, zLanding, ci = "z-transform", 
                 conf.level = 0.95, na.rm = TRUE)
lin.ccc

