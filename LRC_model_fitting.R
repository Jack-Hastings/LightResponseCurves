# Fitting light response curve #
# Harvard Forest 2021 Field collection
# Author: Jack Hastings
# Modified 11/22/2023
# --------------------
# ## Notes ##
# Equations from:
# Marino et al. 2010: https://doi.org/10.1111/j.1365-2435.2009.01630.x
# Schmiege et al. 2022:  https://doi.org/10.1111/pce.14448

# Dependencies:
library(tidyverse)
library(ggplot2)
library(minpack.lm)

# read in light response data
LRC_dataset <- read.csv("HF2021_LRCdata.csv")

# Define Light Response Fit Equations #

# Mitscherlich Equation:
# Q = photosynthetically active photon flux density
# Amax = maximum net assimilation rate
# LCP = light compensation point
# qLCP = apparent quantum yield

Mitscherlich_eqn <- function(Q, Amax, LCP, qLCP){
  
  Amax * (1 - exp((-qLCP * (Q - LCP)) / Amax))
  
}

# Michaelis-Menten equation:
# Gmax = gross maximum photosynthetic rate
# Q = photosynthetically active photon flux density
# K = light intensity when Gmax is half max (half sat)
# Rd = leaf respiration rate

Michaelis_Menten_eqn <- function(Gmax, Q, K, Rd){
  
  ((Gmax * Q)/(K + Q)) - Rd
  
}

# Non-Rectangular Hyperbola equation:
# Q = photosynthetically active photon flux density
# Amax = maximum net assimilation rate
# Rd = leaf respiration rate (dark)
# qLCP = apparent quantum yield
# cF = curvature factor

Non_Rectangular_eqn <- function(qLCP, Q, Amax, cF, Rd){
  
  ((qLCP * Q + Amax - sqrt((qLCP * Q + Amax)**2 -4 * cF * qLCP *Q * Amax)) / 
     (2 * cF)) - Rd

}

# required initial parameter guess for nonlinear least squares models
Mitsch_guess <- c(Amax = 9, LCP = 45, qLCP = 0.05)
Michael_guess <- c(Gmax = 12, K = 175, Rd = 0.5)
NonRect_guess <- c(Amax = 9, qLCP = 0.05, Rd = 0.5, cF = 0.5)



LRC_combo <- LRC_dataset %>% distinct(Species, Tree_Rep)
# Unfortunately, need to manually iterate to ensure convergence
i = 51
Spp <- LRC_combo$Species[i]
Rep <- LRC_combo$Tree_Rep[i]
LRC_ID <- paste0(Spp, "_", Rep)
  
LRC <- LRC_dataset %>%
  filter(Species == Spp & Tree_Rep == Rep)
  

fit_NonRect <- nlsLM(
  formula = A ~ Non_Rectangular_eqn(qLCP, Q, Amax, cF, Rd),
  data = LRC, 
  start = c(qLCP = 0.05, Amax = 9, cF = 0.5, Rd = 0.5), # change this to guess if loop doesn't work
  #lower = c(qLCP = 0.001, Amax = 0, cF = -1.5, Rd = -5),
  #upper = c(qLCP = 0.15, Amax = 50, cF = 1, Rd = 5),
  control = nls.lm.control(maxiter = 1000)

)


fit_NonRect %>% summary %>% coef()

sum_NonRect <- fit_NonRect %>% summary() %>% coef() %>% t() %>%
  as.data.frame() %>% rename(
    Amax_NonRect = Amax,
    qLCP_NonRect = qLCP, 
    Rd_NonRect = Rd, 
    cF_NonRect = cF
  )







