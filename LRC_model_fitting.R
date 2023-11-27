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


LRC_combo <- LRC_dataset %>% distinct(Species, Tree_Rep)
###########
# Need to manually iterate to ensure convergence
###########
i = 19
Spp <- LRC_combo$Species[i]
Rep <- LRC_combo$Tree_Rep[i]
lrc_ID <- paste0(Spp, "_", Rep)
  
LRC <- LRC_dataset %>%
  filter(Species == Spp & Tree_Rep == Rep)
  
# fit nonlinear least square models #

# Mitscherlich:
fit_Mitsch <- nlsLM(
  formula = A ~ Mitscherlich_eqn(Q, Amax, LCP, qLCP),
  data = LRC, 
  start = c(Amax = 9, LCP = 45, qLCP = 0.05),
  control = nls.lm.control(maxiter = 1000)
)

# Michaelis-Menten:
fit_Michael <- nlsLM(
  formula = A ~ Michaelis_Menten_eqn(Gmax, Q, K, Rd),
  data = LRC, 
  start = c(Gmax = 12, K = 175, Rd = 0.5),
  control = nls.lm.control(maxiter = 1000)
)



fit_NonRect <- nlsLM(
  formula = A ~ Non_Rectangular_eqn(qLCP, Q, Amax, cF, Rd),
  data = LRC, 
  start = c(qLCP = 0.05, Amax = 8.6, cF = 0.5, Rd = 2.3), 
  #lower = c(qLCP = 0.001, Amax = 0, cF = -1.5, Rd = -5),
  #upper = c(qLCP = 0.15, Amax = 50, cF = 1.0, Rd = 5),
  control = nls.lm.control(maxiter = 1000)

)

# if NonRect convergence fails, 
# use values from other models as initial parameters
summary(fit_NonRect)
#summary(fit_Michael)
#summary(fit_Mitsch)

# compile the model results

sum_Mitsch <- fit_Mitsch %>% summary %>% coef() %>% t() %>%
  as.data.frame() %>% rename(
    Amax_Mitsch = Amax, 
    LCP_Mitsch = LCP,
    qLCP_Mitsch = qLCP
  )

sum_Michael <- fit_Michael %>% summary %>% coef() %>% t() %>%
  as.data.frame() %>% rename(
    Gmax_Michael = Gmax, 
    K_Micahel = K, 
    Rd_Michael = Rd
  )

sum_NonRect <- fit_NonRect %>% summary() %>% coef() %>% t() %>%
  as.data.frame() %>% rename(
    Amax_NonRect = Amax,
    qLCP_NonRect = qLCP, 
    Rd_NonRect = Rd, 
    cF_NonRect = cF
  )

merged_summary <- bind_cols(sum_Michael, sum_Mitsch, sum_NonRect) %>%
  filter(rownames(.) == 'Estimate') %>%
  mutate(LRC_ID = lrc_ID) %>%
  select(LRC_ID, everything())

#results <- merged_summary # instantiate first time
results <- bind_rows(results, merged_summary)



saveRDS(fit_Mitsch,
        paste0("saved_nls_models/", lrc_ID, "_Mitscherlich.rds"))
saveRDS(fit_Michael,
        paste0("saved_nls_models/", lrc_ID, "_MichaelisMenton.rds"))
saveRDS(fit_NonRect,
        paste0("saved_nls_models/", lrc_ID, "_NonRectHyperbola.rds"))


### Everything below is graphing for visualization

#create a sequence of Q values for plotting
Q_seq <- seq(min(LRC$Q), max(LRC$Q), length.out = 100)

# predict A values using fitted model
predicted_Mitsch <- predict(fit_Mitsch, newdata = data.frame(Q = Q_seq))
predicted_Michael <- predict(fit_Michael, newdata = data.frame(Q = Q_seq))
predicted_NonRect <- predict(fit_NonRect, newdata = data.frame(Q = Q_seq))

# Combine predicted data for plotting
plot_data <- data.frame(Q = Q_seq, 
                        Mitscherlich = predicted_Mitsch,
                        Michaelis_Menten = predicted_Michael,
                        Non_Rectangular = predicted_NonRect)

ggplot(plot_data, aes(x = Q, y = A)) +
  geom_point(data = LRC, aes(color = "Collected Data"),  size = 3) +  
  geom_line(aes(y = Mitscherlich, color = "Mitscherlich"), linetype = "dashed") +  
  geom_line(aes(y = Michaelis_Menten, color = "Michaelis_Menten"), linetype = "dashed") +  
  geom_line(aes(y = Non_Rectangular, color = "Non_Rectangular"), linetype = "dotted") +
  labs(x = "Light Intensity (Q)", y = "Photosynthesis Rate (A)") +
  ggtitle(paste("Light Response Curve -", lrc_ID)) +
  scale_color_manual(
    values = c("Collected Data" = "black", 
               Mitscherlich = "red",
               Michaelis_Menten = "blue",
               Non_Rectangular = "darkgreen")
  ) +
  theme_minimal() +
  theme(
    legend.position = "bottom",
    legend.title = element_blank(),
    legend.spacing.x = unit(0.5, "lines")
  ) +
  guides(
    color = guide_legend(
      title = "Fits",
      title.position = "top",
      title.hjust = 0.5
    )
  )

rm(fit_Michael, fit_Mitsch, fit_NonRect, LRC, merged_summary, sum_Michael,
   sum_Mitsch, sum_NonRect, predicted_Michael, predicted_Mitsch,
   predicted_NonRect, plot_data); gc()

#write.csv(results, "saved_nls_models/model_parameter_summary2.csv")
