################################################################################
### This file is loaded in the `run.R` script to set up the simulation design ##
################################################################################
N              = 1000
# Generate 28 equally spaced follow-up visits ranging from t=0 to 27
mi             = 28
t              = rep(seq(0, mi-1, length.out = mi), N)
logt           = log(t+1, base = mi)
# Generate binary treatment indicator X (`tx`)
pX             = 0.5
# Generate binary baseline status Z
pZ             = 0.5
# Generate 4-level ordinal Y
alpha          = c(-4.75, -1.5, 2)
beta           = c(0, 2, -0.35, 2)  # (X, time, X*time, Z)
gamma.mat      = list()
# gamma setup in the presence of absorbing state (1)
gamma.mat[[1]] = rbind(c( 25, 2, 1),
                       c(-25, 4, 2),
                       c(-25, 1, 4))
# Relax the PO assumption on T only
beta.ppo       = c(-1.75, -4)
ppo.k          = c(2, 3)
