#################################################################
#
#           AUDACE - Summer mortality prediction
#               Part 0: Packages & parameters
#
#################################################################

#------------------------
# Packages
#------------------------

#----- Data management
library(sf) # For cities boundaries
library(tsModel) # For lagging function
library(abind) # To bind array

#----- Analysis
library(splines) # For ns
library(doParallel) # To Paralllelize analysis
library(dlnm) # Perform DLNM
library(FDboost) # For functional regression
library(MASS) # To simulate coefficients

#----- Plotting
library(colorspace)
library(ggplot2)
library(ggmap)
library(rnaturalearthdata)
library(sf)
library(patchwork)
library(fda) # For axis.intervals function

#-------------------------
#   Parameters
#-------------------------

#----- Data

# Data loading
hot_months <- 5:9
clag <- 0:16
climate_labs <- c("AMO", "AO", "NAO", "ONI", "PDO", "PNA", "SOI")
ncl <- length(climate_labs)
region <- c("CMM", "CMQ")
reg_lab <- c("Montreal", "Quebec")

# Colors associated to each climate index
pal <- 1:7

# Labels
longseq <- seq.POSIXt(ISOdate(2000, 1, 1), 
  ISOdate(2002, min(hot_months) - 1, 1), by = "month")
longseq <- longseq[length(longseq) - rev(clag) + 1]
mon_vec <- month.abb[(min(hot_months) - clag - 1) %% 12 + 1]
year_vec <- -((min(hot_months) - clag - 1) %/% 12)

#----- First-stage parameters

# Temp basis parameters
varper <- c(.5, .9)
varfun <- "bs"
vardeg <- 2

# Lag basis parameters
maxlag <- 10
lagnk <- 2

# Degrees of freedom for seasonality and long-term trend
dfseas <- 4
dfdec <- 1

# Heat wave percentile for AN
heatper <- c(.95, .975, .99)

# Persentile for prediction
predper <- c(seq(.1, .9, by = .1), 1:99, seq(99.1, 99.9, by = .1))

# Acceptable percentiles for MMT
mmtper <- c(.1, .9)

# Simulation number for AF confidence intervals
ns <- 1000

#----- Second-stage parameters

# FDboost parameters
nu <- .01
maxstop <- 500

# Confidence intervals
B <- 2000
ci_level <- .95