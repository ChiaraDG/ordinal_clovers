# Longitudinal and ordinal data analyses of cross-sectional, binary outcomes 

The R scripts in this repository contain the code necessary to reproduce simulations and data analysis of the manuscript "Longitudinal and ordinal data analyses of cross-sectional, binary outcomes" (Di Gravio, Tao, Nam, ??, Schildcrout). The repository contains two folder:

* [Simulation Studies](https://github.com/ChiaraDG/ordinal_clovers/tree/main/Simulation%20Studies): contains the code necessary to replicate the simulations in the paper

* [CLOVERS](https://github.com/ChiaraDG/ordinal_clovers/tree/main/CLOVERS): contains the code necessary to replicate the analysis of the CLOVERS study

More detailed instruction on how to reproduce the results are provided in each folder.

## Example

For each simulation scenario, there is a separate R file (`setup.R`) that specifies the parameters used in the data generation. For instance, when generating the data under the OMTM1 scenario:

```
N              = 1000
# 28 equally spaced follow-up visits ranging from t=0 to 27
mi             = 28
t              = rep(seq(0, mi - 1, length.out = mi), N)
logt           = log(t + 1, base = mi)
# binary treatment indicator X (`tx`) and binary baseline status Z
pX             = 0.5
pZ             = 0.5
# 4-level ordinal Y
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
```

The other file (`run.R`) executes a single simulation. It starts with the complete follow-up case, follwed by DAR (mild) and then DAR (severe). For example, when generating the complete follow-up data under the OMTM1 scenario:

```
##### OMTM1 Data - Complete ####
# Generate binary treatment indicator X (`tx`) and binary baseline status Z
tx <- rep(rbinom(N, size = 1, prob = pX), each = mi)
Z <- rep(0, N)
Z[sample(1:N, N * pZ)] <- 1
Z <- rep(Z, each = mi)
# Construct the design matrix `XMat` in the form of (tx, logt, tx*logt, Z)
XMat <- cbind(tx, logt, tx*logt, Z)
# Generate outcome `Y` by relaxing the PO assumption on T
XMat.ppo <- XMat[, c(2, 2)]
id <- rep(1:N, each = mi)
Y <- GenDatOMTM1.ppo(
  id = id, XMat = XMat, alpha = alpha, beta = beta, gamma.mat.list = gamma.mat,
  ppo.k = ppo.k, XMat.ppo = XMat.ppo, beta.ppo = beta.ppo)
# Incorporate the generated data into a dataframe
dat <- data.frame(id = id, X = tx, time = t, logt = logt, Z = Z, Y = Y)
```

Now, the user can decide which modelling procedure to use:

* **Cross-sectional Binary Outcome**

```
# Create binary outcome variable: Ybin=1 (event Y=1), Ybin=0 (no event Y=2:4)
dat$Ybin <- ifelse(dat$Y == 1, 1, 0)
# Subset data at the end of follow-up time only
dat.cross <- subset(dat, logt == 1)

# Fit GLMB model
mdl.glmb <- glm(Ybin ~ X + Z, data = dat.cross, family = "binomial")
```

* **Longitudinal Binary Outcome**

```
# Create binary outcome variable: Ybin=1 (event Y=1), Ybin=2 (no event Y=2:4)
dat$Ybin <- ifelse(dat$Y != 1, 2, 1)

# Fit BMTM model
XMat.dat <- with(dat, cbind(X, logt, X*logt, Z))
Ybin.dat <- dat$Ybin
id.dat <- dat$id
params.bmtm <- c(alpha[1], beta, gamma.mat[[1]][1,1])
mdl.bmtm <- omtm1(params = params.bmtm, yval = Ybin.dat, XMat = XMat.dat, id = id.dat,
                  ppo.k = NULL, x.ppo = NULL, u = NULL, ref.muc = NA,
                  TransIndMtx = NA, ProfileCol = NA, print.level = 1,
                  stepmax = 1, UseGrad = TRUE, iterlim = 250, check.analyticals = FALSE)

# Fit BMTM-C model
XMat.dat.noX <- with(dat, cbind(logt, X*logt, Z))
params.bmtm.noX <- c(alpha[1], beta[-1], gamma.mat[[1]][1, 1])
mdl.bmtm.noX <- omtm1(params = params.bmtm.noX, yval = Ybin.dat, XMat = XMat.dat.noX, id = id.dat,
                      ppo.k = NULL, x.ppo = NULL, u = NULL, ref.muc = NA,
                      TransIndMtx = NA, ProfileCol = NA, print.level = 1,
                      stepmax = 1, UseGrad = TRUE, iterlim = 250, check.analyticals = FALSE)
```

* **Cross-sectional Ordinal Outcome**

```
# Fit GLMO1 model
mdl.glmo1 <- vglm(ordered(Y, levels = 1:4) ~ X + Z, data = dat.cross,
                  family = cumulative(parallel = TRUE))

# Fit GLMO2 model
mat <- cbind(rbind(1, 0, 0), rbind(0, 1, 1))
clist <- list("(Intercept)" = diag(3), "X" = mat, "Z" = rbind(1, 1, 1))
mdl.glmo2 <- vglm(ordered(Y, levels = 1:4) ~ X + Z, data = dat.cross,
                  constraints = clist, family = cumulative(parallel = FALSE))

# Fit GLMO3 model
mdl.glmo3 <- vglm(ordered(Y, levels = 1:4) ~ X + Z, data = dat.cross,
                  family = cumulative(parallel = FALSE ~ X))
```

* **Longitudinal Ordinal Outcome**

```
# Fit OMTM1 model
XMat.ppo1 <- XMat.dat[, c(2, 2)]
Y.dat <- dat$Y
params.omtm1 <- c(alpha, beta, beta.ppo, c(t(gamma.mat[[1]])))
mdl.omtm1 <- omtm1(params = params.omtm1, yval = Y.dat, XMat = XMat.dat, id = id.dat,
                   ppo.k = c(2, 3), x.ppo = XMat.ppo1, u = NULL, ref.muc = NA,
                   TransIndMtx = NA, ProfileCol = NA, print.level = 1,
                   stepmax = 1, UseGrad = TRUE, iterlim = 250, check.analyticals = FALSE)

# Fit OMTM1-C model
XMat.ppo1.noX <- XMat.dat.noX[, c(1, 1)]
params.omtm1.noX <- c(alpha, beta[-1], beta.ppo, c(t(gamma.mat[[1]])))
mdl.omtm1.noX <- omtm1(params = params.omtm1.noX, yval = Y.dat, XMat = XMat.dat.noX, id = id.dat,
                       ppo.k = c(2, 3), x.ppo = XMat.ppo1.noX, u = NULL, ref.muc = NA,
                       TransIndMtx = NA, ProfileCol = NA, print.level = 1,
                       stepmax = 1, UseGrad = TRUE, iterlim = 250, check.analyticals = FALSE)

# Fit OMTM2 model
XMat.ppo2 <- XMat.dat[, c(2, 2, 3)]
Yomtm2.dat <- abs(Y.dat - 5)  # Flip the outcome to fit OMTM2
params.omtm2 <- c(-rev(alpha), -beta[1], -(beta[2] + beta.ppo[2]), -beta[3], -beta[4],
                  beta.ppo[2] - beta.ppo[1], beta.ppo[2], 0,
                  rev(c(t(gamma.mat[[1]]))))
mdl.omtm2 <- omtm1(params = params.omtm2, yval = Yomtm2.dat, XMat = XMat.dat, id = id.dat,
                   ppo.k = c(2, 3, 3), x.ppo = XMat.ppo2, u = NULL, ref.muc = 2,
                   TransIndMtx = NA, ProfileCol = NA, print.level = 1,
                   stepmax = 1, UseGrad = TRUE, iterlim = 250, check.analyticals = FALSE)

# Fit OMTM2-C model
XMat.ppo2.noX <- XMat.dat.noX[, c(1, 1, 2)]
params.omtm2.noX <- c(-rev(alpha), -(beta[2] + beta.ppo[2]), -beta[3], -beta[4],
                      beta.ppo[2] - beta.ppo[1], beta.ppo[2], 0,
                      rev(c(t(gamma.mat[[1]]))))
mdl.omtm2.noX <- omtm1(params = params.omtm2.noX, yval = Yomtm2.dat, XMat = XMat.dat.noX, id = id.dat,
                       ppo.k = c(2, 3, 3), x.ppo = XMat.ppo2.noX, u = NULL, ref.muc = 2,
                       TransIndMtx = NA, ProfileCol = NA, print.level = 1,
                       stepmax = 1, UseGrad = TRUE, iterlim = 250, check.analyticals = FALSE)

# Fit OMTM3 model
XMat.ppo3 <- XMat.dat[, c(2, 2, 3, 3)]
params.omtm3 <- c(alpha, beta, beta.ppo, 0, 0, c(t(gamma.mat[[1]])))
mdl.omtm3 <- omtm1(params = params.omtm3, yval = Y.dat, XMat = XMat.dat, id = id.dat,
                   ppo.k = c(2, 3, 2, 3), x.ppo = XMat.ppo3, u = NULL, ref.muc = NA,
                   TransIndMtx = NA, ProfileCol = NA, print.level = 1,
                   stepmax = 1, UseGrad = TRUE, iterlim = 250, check.analyticals = FALSE)

# Fit OMTM3-C model
XMat.ppo3.noX <- XMat.dat.noX[, c(1, 1, 2, 2)]
params.omtm3.noX <- c(alpha, beta[-1], beta.ppo, 0, 0, c(t(gamma.mat[[1]])))
mdl.omtm3.noX <- omtm1(params = params.omtm3.noX, yval = Y.dat, XMat = XMat.dat.noX, id = id.dat,
                       ppo.k = c(2, 3, 2, 3), x.ppo = XMat.ppo3.noX, u = NULL, ref.muc = NA,
                       TransIndMtx = NA, ProfileCol = NA, print.level = 1,
                       stepmax = stepmax, UseGrad = TRUE, iterlim = 250, check.analyticals = FALSE)
```

* **Survival Outcome**

```
# Create survival data
id.event <- dat[with(dat, logt == 1 & Ybin == 1), ]$id
dat.noevent <- dat %>% filter(!(id %in% id.event)) %>%
  group_by(id) %>% filter(time == max(time)) %>% ungroup()
dat.event <- dat %>% filter(id %in% id.event & Ybin == 1) %>%
  group_by(id) %>% filter(time == min(time)) %>% ungroup()
dat.cox <- rbind(dat.noevent, dat.event)

# Fit Cox PH model
mdl.cox <- coxph(Surv(time, Ybin) ~ X + Z, data = dat.cox, ties = "breslow")
```
