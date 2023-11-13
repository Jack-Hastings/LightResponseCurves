## Fitting light response curves ##
# Harvard Forest 2021
# See Marino et al. 2010 :https://doi.org/10.1111/j.1365-2435.2009.01630.x
# and Schmiege et al. 2022:  https://doi.org/10.1111/pce.14448

# read Schmiege closely -- might be some additional parameters to calculate

library(tidyverse)
library(ggplot2)
library(minpack.lm)

# read in light response data 
LRC_dataset <- read.csv("HF2021_LRCdata.csv")





### Define Mitscherlich equation:
# Q = photosynthetically active photon flux density
# Amax = maximum net assimilation rate 
# LCP = light compensation point
# qLCP = apparent quantum yield

Mitscherlich_eqn <- function(Q, Amax, LCP, qLCP){
  
  Amax * (1 - exp((-qLCP * (Q - LCP))/Amax))
  
}


### Define Michaelis-Menten equation:
# Gmax = gross maximum photosynthetic rate
# Q = photosynthetically active photon flux density
# K = light intensity when Gmax is half max (half sat)
# Rd = leaf respiration rate

Michaelis_Menten_eqn <- function(Gmax, Q, K, Rd){
  
  ((Gmax * Q)/(K + Q)) - Rd
  
}

### Define Non-Rectangular Hyperbola equation:
# Q = photosynthetically active photon flux density
# Amax = maximum net assimilation rate
# Rd = leaf respiration rate
# qLCP = apparent quantum yield
# cF = curvature factor


Non_Rectangular_eqn <- function(qLCP, Q, Amax, cF, Rd){
  
  ((qLCP*Q + Amax - sqrt((qLCP*Q + Amax)**2 - 4*cF*qLCP*Q*Amax)) /
     (2*cF)) - Rd
}



#initial parameter guess for nls models
Mitsch_guess <-c(Amax = 9, LCP = 45, qLCP = 0.05)
Michael_guess <- c(Gmax = 12, K = 175, Rd = 2) 
NonRect_guess <- c(Amax = 3.38, qLCP = 0.0035, Rd = 0.5, cF = 0.50)




######################################################################
# iterate from here...manually...
# non-rectangular hyperbola is finicky and needs some manual input...
######################################################################



#Pull out individual LRC curve based on species and tree rep ID
unique_combinations <- LRC_dataset %>% distinct(Species, Tree_Rep)

i = 51
#51 failed -- come back to it...can't figure it out. 

Spp <-unique_combinations$Species[i]
Rep <- unique_combinations$Tree_Rep[i]

LRC <- LRC_dataset %>%
  filter(Species == Spp & Tree_Rep == Rep)

lrc_ID <- paste0(Spp, "_", Rep)

plot(LRC$A~LRC$Q)


# Fit the models using nls
#Mitscherlich:
fit_Mitsch <- nls(A ~ Mitscherlich_eqn(Q, Amax, LCP, qLCP), 
                  data = LRC, start = Mitsch_guess)


#Michaelis-Menten:
fit_Michael <- nls(A ~ Michaelis_Menten_eqn(Gmax, Q, K, Rd),
                   data = LRC, start = Michael_guess)

#Non-Rectangular Hyperbola:
# note -- run w/ algorithm = "default". If that fails, run w/ algorithm = "port"
# if that fails, tweak initial guess parameters
fit_NonRect <- nls(A ~ Non_Rectangular_eqn(qLCP, Q, Amax, cF, Rd),
                   data = LRC, start = NonRect_guess, 
                   upper = c(Amax = 50, qLCP = 0.15, Rd = 5, cF = 1),
                   lower = c(Amax = -Inf, qLCP = 0.001, Rd = -5, cF = -1.5)
                   ,algorithm = "port",
                   trace = TRUE,
                   warnOnly = TRUE
                   )





##############
# Define your model function
Non_Rectangular_eqn_2 <- function(qLCP, Q, Amax, cF, Rd) {
  
  ((qLCP * Q + Amax - sqrt((qLCP * Q + Amax)^2 - 4 * cF * qLCP * Q * Amax)) / (2 * cF)) - Rd
}

# Define initial parameter guesses for all parameters
NonRect_guess <- c(qLCP = 0.01, Amax = 10, cF = 0.5, Rd = 1)

# Fit the model using nlsLM from minpack.lm
fit_NonRect <- nlsLM(
  formula = A ~ Non_Rectangular_eqn_2(qLCP, Q, Amax, cF, Rd),
  data = LRC,
  start = NonRect_guess,
  lower = c(qLCP = 0.001, Amax = 0, cF = -1.5, Rd = -5),
  upper = c(qLCP = 0.15, Amax = 50, cF = 1, Rd = 5),
  control = nls.lm.control(maxiter = 1000)  # Adjust maxiter as needed
)

####





sum_Mitsch <- fit_Mitsch %>% summary() %>% coef() %>% t() %>%
  as.data.frame() %>% rename(
    Amax_Mitsch = Amax, 
    LCP_Mitsch = LCP, 
    qLCP_Mitsch = qLCP)

sum_Michael <- fit_Michael %>% summary() %>% coef() %>% t() %>% 
  as.data.frame() %>% rename(
    Gmax_Michael = Gmax, K_Michael = K, Rd_Michael = Rd)

sum_NonRect <- fit_NonRect %>% summary() %>% coef() %>% t() %>% 
  as.data.frame %>% rename(
    Amax_NonRect = Amax, qLCP_NonRect = qLCP, RD_NonRect = Rd, cF_NonRect = cF)
    

merged_summary <- bind_cols(sum_Michael, sum_Mitsch, sum_NonRect) %>%
  filter(rownames(.) == 'Estimate') %>% 
  mutate(LRC_ID = lrc_ID) %>% 
  select(LRC_ID, everything())


#results <- merged_summary #just first time
results <- bind_rows(results, merged_summary)



saveRDS(fit_Mitsch,
        paste0("saved_nls_models/", lrc_ID, "_Mitscherlich.rds"))
saveRDS(fit_Michael,
        paste0("saved_nls_models/", lrc_ID, "_MichaelisMenton.rds"))
saveRDS(fit_NonRect,
        paste0("saved_nls_models/", lrc_ID, "_NonRectHyperbola.rds"))



## fit the model and plot the results

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






