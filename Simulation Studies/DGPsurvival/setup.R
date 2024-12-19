################################################################################
### This file is loaded in the `run.R` script to set up the simulation design ##
################################################################################
N         = 1000
# Generate 28 equally spaced follow-up visits ranging from t=0 to 27
mi        = 28
t         = rep(seq(0, mi-1, length.out = mi), N)
# Generate binary treatment indicator X (`tx`)
pX        = 0.5
# Generate binary baseline status Z 
pZ        = 0.5
# Survival data parameters (Weibull distribution)
lambda    = 0.01
nu        = 0.5
beta      = c(-0.35, 2)  # (X, Z)
