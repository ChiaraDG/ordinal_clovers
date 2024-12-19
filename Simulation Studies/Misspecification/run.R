##### This file executes a single simulation and saves data and results
# Load packages
library(OMTM1)  
library(dplyr)
library(survival)
library(VGAM)
# Check the version of loaded packages
print(packageVersion("OMTM1"));
print(packageVersion("dplyr"));
print(packageVersion("survival"))
print(packageVersion("VGAM"))
# Load simulation design setup
source("setup.R")

simulation <- function(run){
  
  ##### Data - Complete (misspecification of beta_xt) ####
  cat("Generating random data from known parameters\n")
  # Generate binary treatment indicator X (`tx`)
  tx <- rep(rbinom(N, size = 1, prob = pX), each = mi)
  # Generate binary baseline status Z (Fix p(Z) = 0.5)
  Z <- rep(0, N)
  Z[sample(1:N, N * pZ)] <- 1
  Z <- rep(Z, each = mi)
  # Construct the design matrix `XMat` in the form of (tx, logt, tx*logt, Z)
  XMat <- cbind(tx, logt, tx*logt, Z)
  # Generate outcome `Y` by relaxing the PO assumption on T and XT
  # in the presence of the absorbing state (k=1)
  XMat.ppo <- XMat[, c(2, 2, 3, 3)]
  beta.ppo <- c(-1.75, -4, -0.2, -0.2)
  ppo.k <- c(2, 3, 2, 3)
  id <- rep(1:N, each = mi)
  Y <- GenDatOMTM1.ppo(
    id = id, XMat = XMat, alpha = alpha, beta = beta, gamma.mat.list = gamma.mat,
    ppo.k = ppo.k, XMat.ppo = XMat.ppo, beta.ppo = beta.ppo)
  # Incorporate the generated data into a dataframe
  dat <- data.frame(id = id, X = tx, time = t, logt = logt, Z = Z, Y = Y)
  # Create binary outcome variable: Ybin=1 (event Y=1), Ybin=0 (no event Y=2:4)
  dat$Ybin <- ifelse(dat$Y == 1, 1, 0)
  dat0 <- dat
  id.event <- dat[with(dat, logt == 1 & Ybin == 1), ]$id
  dat.noevent <- dat %>% filter(!(id %in% id.event)) %>%
    group_by(id) %>% filter(time == max(time)) %>% ungroup()
  dat.event <- dat %>% filter(id %in% id.event & Ybin == 1) %>%
    group_by(id) %>% filter(time == min(time)) %>% ungroup()
  dat.cox <- rbind(dat.noevent, dat.event)
  # Save data
  if (file.exists("data-missX")){
    save(dat, file = paste0("data-missX/data-", run, ".RData"))}

  # Calculate p(Z) at the beginning of the study
  dat.t0 <- dat %>% filter(logt==0)
  pZt0 <- mean(dat.t0$Z==1)

  ##### GLMB ####
  cat("Fitting GLMB\n")
  # Subset data at the end of follow-up time only
  dat.cross <- subset(dat, logt == 1)
  dat.cross$Ybin <- ifelse(dat.cross$Ybin == 1, 1, 0)
  # Fit GLMB model and extract coefficients, covariance matrix
  mdl.glmb <- glm(Ybin ~ X + Z, data = dat.cross, family = "binomial")
  coefs.glmb <- coef(mdl.glmb)  # (beta_0, beta_x, beta_z)
  names(coefs.glmb) <- c("est.a1", "est.bx1", "est.bz")
  vcov.glmb <- vcov(mdl.glmb)
  ses.glmb <- sqrt(diag(vcov.glmb))
  names(ses.glmb) <- c("se.a1", "se.bx1", "se.bz")
  # Calculate probabilities
  pX1Z1 <- expit(sum(coefs.glmb))
  pX1Z0 <- expit(sum(coefs.glmb[c(1, 2)]))
  pX0Z1 <- expit(sum(coefs.glmb[c(1, 3)]))
  pX0Z0 <- expit(sum(coefs.glmb[1]))
  # Calculate the Jacobian matrix
  J.glmb <- matrix(c(
    ((1 - pX1Z1) * pX1Z1 * pZt0 - (1 - pX0Z1) * pX0Z1 * pZt0) +
      ((1 - pX1Z0) * pX1Z0 * (1 - pZt0) - (1 - pX0Z0) * pX0Z0 * (1 - pZt0)),
    (1 - pX1Z1) * pX1Z1 * pZt0 + (1 - pX1Z0) * pX1Z0 * (1 - pZt0),
    (1 - pX1Z1) * pX1Z1 * pZt0 - (1 - pX0Z1) * pX0Z1 * pZt0), nrow = 1)
  # Calculate RD estimate and its standard error
  RDest.glmb <- (pX1Z1 * pZt0 - pX0Z1 * pZt0) +
    (pX1Z0 * (1 - pZt0) - pX0Z0 * (1 - pZt0))
  RDse.glmb <- c(sqrt(J.glmb %*% vcov.glmb %*% t(J.glmb)))

  ##### BMTM ####
  cat("Fitting BMTM\n")
  # Create binary outcome variable: Ybin=1 (event Y=1), Ybin=2 (no event Y=2:4)
  dat$Ybin <- ifelse(dat$Y != 1, 2, 1)
  # Fit BMTM model and extract coefficients, covariance matrix
  XMat.dat <- with(dat, cbind(X, logt, X*logt, Z))
  Ybin.dat <- dat$Ybin
  id.dat <- dat$id
  params.bmtm <- c(alpha[1], beta, gamma.mat[[1]][1, 1])
  stepmax <- 1
  code <- 5
  while(code == 5){
    cat(paste("Try stepmax =", stepmax, "\n"))
    reset_step <- FALSE
    tryCatch(
      mdl.bmtm <- omtm1(
        params = params.bmtm, yval = Ybin.dat, XMat = XMat.dat, id = id.dat,
        ppo.k = NULL, x.ppo = NULL, u = NULL, ref.muc = NA,
        TransIndMtx = NA, ProfileCol = NA, print.level = 1,
        stepmax = stepmax, UseGrad = TRUE, iterlim = 250, check.analyticals = FALSE),
      error = function(e){reset_step <<- TRUE})
    cat(paste("Reset step", reset_step, "\n"))
    if(reset_step == TRUE){
      stepmax <- stepmax * 0.25
    } else {
      stepmax <- stepmax * 1.5
      params.bmtm <- c(mdl.bmtm$alpha, mdl.bmtm$beta, mdl.bmtm$gamma[[1]][1, 1])
      code <- code.bmtm <- mdl.bmtm$control["convergence_code"]
    }}
  coefs.bmtm <- c(mdl.bmtm$alpha, mdl.bmtm$beta, mdl.bmtm$gamma[[1]][1, 1])
  names(coefs.bmtm) <- c("est.a1", "est.bx1", "est.bt1", "est.bxt1", "est.bz", "est.g11")
  vcov.bmtm <- mdl.bmtm$vcov
  ses.bmtm <- c(sqrt(diag(vcov.bmtm)), vcov.bmtm[2, 4])
  names(ses.bmtm) <- c("se.a1", "se.bx1", "se.bt1", "se.bxt1", "se.bz", "cov.bx1bxt1")
  # Calculate probabilities
  pX1Z1 <- expit(sum(coefs.bmtm[paste0("est.", c("a1", "bx1", "bt1", "bxt1", "bz"))]))
  pX1Z0 <- expit(sum(coefs.bmtm[paste0("est.", c("a1", "bx1", "bt1", "bxt1"))]))
  pX0Z1 <- expit(sum(coefs.bmtm[paste0("est.", c("a1", "bt1", "bz"))]))
  pX0Z0 <- expit(sum(coefs.bmtm[paste0("est.", c("a1", "bt1"))]))
  # Calculate the Jacobian matrix
  J0 <- Jt <- ((1 - pX1Z1) * pX1Z1 * pZt0 - (1 - pX0Z1) * pX0Z1 * pZt0) +
    ((1 - pX1Z0) * pX1Z0 * (1 - pZt0) - (1 - pX0Z0) * pX0Z0 * (1 - pZt0))
  Jx <- Jxt <- (1 - pX1Z1) * pX1Z1 * pZt0 + (1 - pX1Z0) * pX1Z0 * (1 - pZt0)
  Jz <- (1 - pX1Z1) * pX1Z1 * pZt0 - (1 - pX0Z1) * pX0Z1 * pZt0
  J.bmtm <- matrix(c(J0, Jx, Jt, Jxt, Jz), nrow = 1)
  RDest.bmtm <- (pX1Z1 * pZt0 - pX0Z1 * pZt0) +
    (pX1Z0 * (1 - pZt0) - pX0Z0 * (1 - pZt0))
  RDse.bmtm <- c(sqrt(J.bmtm %*% vcov.bmtm %*% t(J.bmtm)))
  cat(paste("BMTM convergence code =", code.bmtm, "\n"))

  ##### BMTM-C ####
  cat("Fitting BMTM-C\n")
  # Fit BMTM-C model and extract coefficients, covariance matrix
  XMat.dat.noX <- with(dat, cbind(logt, X*logt, Z))
  params.bmtm.noX <- c(alpha[1], beta[-1], gamma.mat[[1]][1, 1])
  stepmax <- 1
  code <- 5
  while(code == 5){
    cat(paste("Try stepmax =", stepmax, "\n"))
    reset_step <- FALSE
    tryCatch(
      mdl.bmtm.noX <- omtm1(
        params = params.bmtm.noX, yval = Ybin.dat, XMat = XMat.dat.noX, id = id.dat,
        ppo.k = NULL, x.ppo = NULL, u = NULL, ref.muc = NA,
        TransIndMtx = NA, ProfileCol = NA, print.level = 1,
        stepmax = 1, UseGrad = TRUE, iterlim = 250, check.analyticals = FALSE),
      error = function(e){reset_step <<- TRUE})
    cat(paste("Reset step", reset_step, "\n"))
    if(reset_step == TRUE){
      stepmax <- stepmax * 0.25
    } else {
      stepmax <- stepmax * 1.5
      params.bmtm.noX <- c(mdl.bmtm.noX$alpha, mdl.bmtm.noX$beta, mdl.bmtm.noX$gamma[[1]][1, 1])
      code <- code.bmtm.noX <- mdl.bmtm.noX$control["convergence_code"]
    }}
  coefs.bmtm.noX <- c(mdl.bmtm.noX$alpha, mdl.bmtm.noX$beta, mdl.bmtm.noX$gamma[[1]][1, 1])
  names(coefs.bmtm.noX) <- c("est.a1", "est.bt1", "est.bxt1", "est.bz", "est.g11")
  vcov.bmtm.noX <- mdl.bmtm.noX$vcov
  ses.bmtm.noX <- sqrt(diag(vcov.bmtm.noX))
  names(ses.bmtm.noX) <- c("se.a1", "se.bt1", "se.bxt1", "se.bz")
  # Calculate probabilities
  pX1Z1 <- expit(sum(coefs.bmtm.noX[paste0("est.", c("a1", "bt1", "bxt1", "bz"))]))
  pX1Z0 <- expit(sum(coefs.bmtm.noX[paste0("est.", c("a1", "bt1", "bxt1"))]))
  pX0Z1 <- expit(sum(coefs.bmtm.noX[paste0("est.", c("a1", "bt1", "bz"))]))
  pX0Z0 <- expit(sum(coefs.bmtm.noX[paste0("est.", c("a1", "bt1"))]))
  # Calculate the Jacobian matrix
  J0 <- Jt <- ((1 - pX1Z1) * pX1Z1 * pZt0 - (1 - pX0Z1) * pX0Z1 * pZt0) +
    ((1 - pX1Z0) * pX1Z0 * (1 - pZt0) - (1 - pX0Z0) * pX0Z0 * (1 - pZt0))
  Jxt <- (1 - pX1Z1) * pX1Z1 * pZt0 + (1 - pX1Z0) * pX1Z0 * (1 - pZt0)
  Jz <- (1 - pX1Z1) * pX1Z1 * pZt0 - (1 - pX0Z1) * pX0Z1 * pZt0
  J.bmtm.noX <- matrix(c(J0, Jt, Jxt, Jz), nrow = 1)
  RDest.bmtm.noX <- (pX1Z1 * pZt0 - pX0Z1 * pZt0) +
    (pX1Z0 * (1 - pZt0) - pX0Z0 * (1 - pZt0))
  RDse.bmtm.noX <- c(sqrt(J.bmtm.noX %*% vcov.bmtm.noX %*% t(J.bmtm.noX)))
  cat(paste("BMTM-C convergence code =", code.bmtm.noX, "\n"))

  ##### GLMO1 #####
  # Fit GLMB model and extract coefficients, covariance matrix
  mdl.glmo1 <- vglm(ordered(Y, levels = 1:4) ~ X + Z, data = dat.cross,
                    family = cumulative(parallel = TRUE))
  coefs.glmo1 <- coef(mdl.glmo1)
  names(coefs.glmo1) <- c("est.a1", "est.a2", "est.a3", "est.bx1", "est.bz")
  vcov.glmo1 <- vcov(mdl.glmo1)
  ses.glmo1 <- sqrt(diag(vcov.glmo1))
  names(ses.glmo1) <- c("se.a1", "se.a2", "se.a3", "se.bx1", "se.bz")
  vcov.glmo1 <- vcov.glmo1[c(1, 4:5),c(1, 4:5)]
  # Calculate probabilities
  pX1Z1 <- expit(sum(coefs.glmo1[paste0("est.", c("a1", "bx1", "bz"))]))
  pX1Z0 <- expit(sum(coefs.glmo1[paste0("est.", c("a1", "bx1"))]))
  pX0Z1 <- expit(sum(coefs.glmo1[paste0("est.", c("a1", "bz"))]))
  pX0Z0 <- expit(sum(coefs.glmo1[paste0("est.", c("a1"))]))
  # Calculate the Jacobian matrix
  J.glmo1 <- matrix(c(
    ((1 - pX1Z1) * pX1Z1 * pZt0 - (1 - pX0Z1) * pX0Z1 * pZt0) +
      ((1 - pX1Z0) * pX1Z0 * (1 - pZt0) - (1 - pX0Z0) * pX0Z0 * (1 - pZt0)),
    (1 - pX1Z1) * pX1Z1 * pZt0 + (1 - pX1Z0) * pX1Z0 * (1 - pZt0),
    (1 - pX1Z1) * pX1Z1 * pZt0 - (1 - pX0Z1) * pX0Z1 * pZt0), nrow = 1)
  # Calculate RD estimate and its standard error
  RDest.glmo1 <- (pX1Z1 * pZt0 - pX0Z1 * pZt0) +
    (pX1Z0 * (1 - pZt0) - pX0Z0 * (1 - pZt0))
  RDse.glmo1 <- c(sqrt(J.glmo1 %*% vcov.glmo1 %*% t(J.glmo1)))

  ##### GLMO2 #####
  # Fit GLMB model and extract coefficients, covariance matrix
  mat <- cbind(rbind(1, 0, 0), rbind(0, 1, 1))
  clist <- list("(Intercept)" = diag(3), "X" = mat, "Z" = rbind(1, 1, 1))
  mdl.glmo2 <- vglm(ordered(Y, levels = 1:4) ~ X + Z, data = dat.cross,
                    constraints = clist, family = cumulative(parallel = FALSE))
  coefs.glmo2 <- coef(mdl.glmo2)
  names(coefs.glmo2) <- c("est.a1", "est.a2", "est.a3", "est.bx1", "est.bx2", "est.bz")
  vcov.glmo2 <- vcov(mdl.glmo2)
  ses.glmo2 <- sqrt(diag(vcov.glmo2))
  names(ses.glmo2) <- c("se.a1", "se.a2", "se.a3", "se.bx1", "se.bx2", "se.bz")
  vcov.glmo2 <- vcov.glmo2[c(1, 4, 6), c(1, 4, 6)]
  # Calculate probabilities
  pX1Z1 <- expit(sum(coefs.glmo2[paste0("est.", c("a1", "bx1", "bz"))]))
  pX1Z0 <- expit(sum(coefs.glmo2[paste0("est.", c("a1", "bx1"))]))
  pX0Z1 <- expit(sum(coefs.glmo2[paste0("est.", c("a1", "bz"))]))
  pX0Z0 <- expit(sum(coefs.glmo2[paste0("est.", c("a1"))]))
  # Calculate the Jacobian matrix
  J.glmo2 <- matrix(c(
    ((1 - pX1Z1) * pX1Z1 * pZt0 - (1 - pX0Z1) * pX0Z1 * pZt0) +
      ((1 - pX1Z0) * pX1Z0 * (1 - pZt0) - (1 - pX0Z0) * pX0Z0 * (1 - pZt0)),
    (1 - pX1Z1) * pX1Z1 * pZt0 + (1 - pX1Z0) * pX1Z0 * (1 - pZt0),
    (1 - pX1Z1) * pX1Z1 * pZt0 - (1 - pX0Z1) * pX0Z1 * pZt0), nrow = 1)
  # Calculate RD estimate and its standard error
  RDest.glmo2 <- (pX1Z1 * pZt0 - pX0Z1 * pZt0) +
    (pX1Z0 * (1 - pZt0) - pX0Z0 * (1 - pZt0))
  RDse.glmo2 <- c(sqrt(J.glmo2 %*% vcov.glmo2 %*% t(J.glmo2)))

  ##### GLMO3 #####
  # Fit GLMB model and extract coefficients, covariance matrix
  mdl.glmo3 <- vglm(ordered(Y, levels = 1:4) ~ X + Z, data = dat.cross,
                    family = cumulative(parallel = FALSE ~ X))
  coefs.glmo3 <- coef(mdl.glmo3)
  names(coefs.glmo3) <- c("est.a1", "est.a2", "est.a3", "est.bx1", "est.bx2", "est.bx3", "est.bz")
  vcov.glmo3 <- vcov(mdl.glmo3)
  ses.glmo3 <- sqrt(diag(vcov.glmo3))
  names(ses.glmo3) <- c("se.a1", "se.a2", "se.a3", "se.bx1", "se.bx2", "se.bx3", "se.bz")
  vcov.glmo3 <- vcov.glmo3[c(1, 4, 7), c(1, 4, 7)]
  # Calculate probabilities
  pX1Z1 <- expit(sum(coefs.glmo3[paste0("est.", c("a1", "bx1", "bz"))]))
  pX1Z0 <- expit(sum(coefs.glmo3[paste0("est.", c("a1", "bx1"))]))
  pX0Z1 <- expit(sum(coefs.glmo3[paste0("est.", c("a1", "bz"))]))
  pX0Z0 <- expit(sum(coefs.glmo3[paste0("est.", c("a1"))]))
  # Calculate the Jacobian matrix
  J.glmo3 <- matrix(c(
    ((1 - pX1Z1) * pX1Z1 * pZt0 - (1 - pX0Z1) * pX0Z1 * pZt0) +
      ((1 - pX1Z0) * pX1Z0 * (1 - pZt0) - (1 - pX0Z0) * pX0Z0 * (1 - pZt0)),
    (1 - pX1Z1) * pX1Z1 * pZt0 + (1 - pX1Z0) * pX1Z0 * (1 - pZt0),
    (1 - pX1Z1) * pX1Z1 * pZt0 - (1 - pX0Z1) * pX0Z1 * pZt0), nrow = 1)
  # Calculate RD estimate and its standard error
  RDest.glmo3 <- (pX1Z1 * pZt0 - pX0Z1 * pZt0) +
    (pX1Z0 * (1 - pZt0) - pX0Z0 * (1 - pZt0))
  RDse.glmo3 <- c(sqrt(J.glmo3 %*% vcov.glmo3 %*% t(J.glmo3)))

  ##### OMTM1 #####
  cat("Fitting OMTM1\n")
  XMat.ppo1 <- XMat.dat[, c(2, 2)]
  Y.dat <- dat$Y
  params.omtm1 <- c(alpha, beta, beta.ppo[1:2], c(t(gamma.mat[[1]])))
  stepmax <- 1
  code <- 5
  while(code == 5){
    cat(paste("Try stepmax =", stepmax, "\n"))
    reset_step <- FALSE
    tryCatch(
      mdl.omtm1 <- omtm1(
        params = params.omtm1, yval = Y.dat, XMat = XMat.dat, id = id.dat,
        ppo.k = c(2, 3), x.ppo = XMat.ppo1, u = NULL, ref.muc = NA,
        TransIndMtx = NA, ProfileCol = NA, print.level = 1,
        stepmax = stepmax, UseGrad = TRUE, iterlim = 250, check.analyticals = FALSE),
      error = function(e){reset_step <<- TRUE})
    cat(paste("Reset step", reset_step, "\n"))
    if(reset_step == TRUE){
      stepmax <- stepmax * 0.25
    } else {
      stepmax <- stepmax * 1.5
      params.omtm1 <- c(mdl.omtm1$alpha, mdl.omtm1$beta, c(t(mdl.omtm1$gamma[[1]])))
      code <- code.omtm1 <- mdl.omtm1$control["convergence_code"]
    }}
  coefs.omtm1 <- c(mdl.omtm1$alpha, mdl.omtm1$beta, c(t(mdl.omtm1$gamma[[1]])))
  names(coefs.omtm1) <- c("est.a1", "est.a2", "est.a3",
                          "est.bx1", "est.bt1", "est.bxt1", "est.bz",
                          "est.bt2", "est.bt3",
                          "est.g11", "est.g12", "est.g13",
                          "est.g21", "est.g22", "est.g23",
                          "est.g31", "est.g32", "est.g33")
  vcov.omtm1 <- mdl.omtm1$vcov
  ses.omtm1 <- c(sqrt(diag(vcov.omtm1)), vcov.omtm1[4, 6])
  names(ses.omtm1) <- c("se.a1", "se.a2", "se.a3",
                        "se.bx1", "se.bt1", "se.bxt1", "se.bz",
                        "se.bt2", "se.bt3",
                        "se.g12", "se.g13",
                        "se.g22", "se.g23",
                        "se.g32", "se.g33", "cov.bx1bxt1")
  vcov.omtm1 <- vcov.omtm1[c(1, 4:7), c(1, 4:7)]
  # Calculate probabilities
  pX1Z1 <- expit(sum(coefs.omtm1[paste0("est.", c("a1", "bx1", "bt1", "bxt1", "bz"))]))
  pX1Z0 <- expit(sum(coefs.omtm1[paste0("est.", c("a1", "bx1", "bt1", "bxt1"))]))
  pX0Z1 <- expit(sum(coefs.omtm1[paste0("est.", c("a1", "bt1", "bz"))]))
  pX0Z0 <- expit(sum(coefs.omtm1[paste0("est.", c("a1", "bt1"))]))
  # Calculate the Jacobian matrix
  J0 <- Jt <- ((1 - pX1Z1) * pX1Z1 * pZt0 - (1 - pX0Z1) * pX0Z1 * pZt0) +
    ((1 - pX1Z0) * pX1Z0 * (1 - pZt0) - (1 - pX0Z0) * pX0Z0 * (1 - pZt0))
  Jx <- Jxt <- (1 - pX1Z1) * pX1Z1 * pZt0 + (1 - pX1Z0) * pX1Z0 * (1 - pZt0)
  Jz <- (1 - pX1Z1) * pX1Z1 * pZt0 - (1 - pX0Z1) * pX0Z1 * pZt0
  J.omtm1 <- matrix(c(J0, Jx, Jt, Jxt, Jz), nrow = 1)
  RDest.omtm1 <- (pX1Z1 * pZt0 - pX0Z1 * pZt0) +
    (pX1Z0 * (1 - pZt0) - pX0Z0 * (1 - pZt0))
  RDse.omtm1 <- c(sqrt(J.omtm1 %*% vcov.omtm1 %*% t(J.omtm1)))
  cat(paste("OMTM1 convergence code =", code.omtm1, "\n"))

  ##### OMTM1-C #####
  cat("Fitting OMTM1-C\n")
  XMat.ppo1.noX <- XMat.dat.noX[, c(1, 1)]
  params.omtm1.noX <- c(alpha, beta[-1], beta.ppo[1:2], c(t(gamma.mat[[1]])))
  stepmax <- 1
  code <- 5
  while(code == 5){
    cat(paste("Try stepmax =", stepmax, "\n"))
    reset_step <- FALSE
    tryCatch(
      mdl.omtm1.noX <- omtm1(
        params = params.omtm1.noX, yval = Y.dat, XMat = XMat.dat.noX, id = id.dat,
        ppo.k = c(2, 3), x.ppo = XMat.ppo1.noX, u = NULL, ref.muc = NA,
        TransIndMtx = NA, ProfileCol = NA, print.level = 1,
        stepmax = stepmax, UseGrad = TRUE, iterlim = 250, check.analyticals = FALSE),
      error = function(e){reset_step <<- TRUE})
    cat(paste("Reset step", reset_step, "\n"))
    if(reset_step == TRUE){
      stepmax <- stepmax * 0.25
    } else {
      stepmax <- stepmax * 1.5
      params.omtm1.noX <- c(mdl.omtm1.noX$alpha, mdl.omtm1.noX$beta, c(t(mdl.omtm1.noX$gamma[[1]])))
      code <- code.omtm1.noX <- mdl.omtm1.noX$control["convergence_code"]
    }}
  coefs.omtm1.noX <- c(mdl.omtm1.noX$alpha, mdl.omtm1.noX$beta, c(t(mdl.omtm1.noX$gamma[[1]])))
  names(coefs.omtm1.noX) <- c("est.a1", "est.a2", "est.a3",
                              "est.bt1", "est.bxt1", "est.bz",
                              "est.bt2", "est.bt3",
                              "est.g11", "est.g12", "est.g13",
                              "est.g21", "est.g22", "est.g23",
                              "est.g31", "est.g32", "est.g33")
  vcov.omtm1.noX <- mdl.omtm1.noX$vcov
  ses.omtm1.noX <- sqrt(diag(vcov.omtm1.noX))
  names(ses.omtm1.noX) <- c("se.a1", "se.a2", "se.a3",
                            "se.bt1", "se.bxt1", "se.bz",
                            "se.bt2", "se.bt3",
                            "se.g12", "se.g13",
                            "se.g22", "se.g23",
                            "se.g32", "se.g33")
  vcov.omtm1.noX <- vcov.omtm1.noX[c(1, 4:6), c(1, 4:6)]
  # Calculate probabilities
  pX1Z1 <- expit(sum(coefs.omtm1.noX[paste0("est.", c("a1", "bt1", "bxt1", "bz"))]))
  pX1Z0 <- expit(sum(coefs.omtm1.noX[paste0("est.", c("a1", "bt1", "bxt1"))]))
  pX0Z1 <- expit(sum(coefs.omtm1.noX[paste0("est.", c("a1", "bt1", "bz"))]))
  pX0Z0 <- expit(sum(coefs.omtm1.noX[paste0("est.", c("a1", "bt1"))]))
  # Calculate the Jacobian matrix
  J0 <- Jt <- ((1 - pX1Z1) * pX1Z1 * pZt0 - (1 - pX0Z1) * pX0Z1 * pZt0) +
    ((1 - pX1Z0) * pX1Z0 * (1 - pZt0) - (1 - pX0Z0) * pX0Z0 * (1 - pZt0))
  Jxt <- (1 - pX1Z1) * pX1Z1 * pZt0 + (1 - pX1Z0) * pX1Z0 * (1 - pZt0)
  Jz <- (1 - pX1Z1) * pX1Z1 * pZt0 - (1 - pX0Z1) * pX0Z1 * pZt0
  J.omtm1.noX <- matrix(c(J0, Jt, Jxt, Jz), nrow = 1)
  RDest.omtm1.noX <- (pX1Z1 * pZt0 - pX0Z1 * pZt0) +
    (pX1Z0 * (1 - pZt0) - pX0Z0 * (1 - pZt0))
  RDse.omtm1.noX <- c(sqrt(J.omtm1.noX %*% vcov.omtm1.noX %*% t(J.omtm1.noX)))
  cat(paste("OMTM1-C convergence code =", code.omtm1.noX, "\n"))

  ##### OMTM2 #####
  cat("Fitting OMTM2\n")
  XMat.ppo2 <- XMat.dat[, c(2, 2, 3)]
  # Flip the outcome to fit OMTM2
  Yomtm2.dat <- abs(Y.dat - 5)
  params.omtm2 <- c(-rev(alpha), -beta[1], -(beta[2] + beta.ppo[2]), -beta[3], -beta[4],
                    beta.ppo[2] - beta.ppo[1], beta.ppo[2], 0,
                    rev(c(t(gamma.mat[[1]]))))
  stepmax <- 1
  code <- 5
  while(code == 5){
    cat(paste("Try stepmax =", stepmax, "\n"))
    reset_step <- FALSE
    tryCatch(mdl.omtm2 <- omtm1(
      params = params.omtm2, yval = Yomtm2.dat, XMat = XMat.dat, id = id.dat,
      ppo.k = c(2, 3, 3), x.ppo = XMat.ppo2, u = NULL, ref.muc = 2,
      TransIndMtx = NA, ProfileCol = NA, print.level = 1,
      stepmax = stepmax, UseGrad = TRUE, iterlim = 250, check.analyticals = FALSE),
      error = function(e){reset_step <<- TRUE})
    cat(paste("Reset step", reset_step, "\n"))
    if(reset_step == TRUE){
      stepmax <- stepmax * 0.25
    } else {
      stepmax <- stepmax * 1.5
      params.omtm2 <- c(mdl.omtm2$alpha, mdl.omtm2$beta, c(t(mdl.omtm2$gamma[[1]])))
      code <- code.omtm2 <- mdl.omtm2$control["convergence_code"]
    }}
  coefs.omtm2 <- c(mdl.omtm2$alpha, mdl.omtm2$beta, c(t(mdl.omtm2$gamma[[1]])))
  names(coefs.omtm2) <- c("est.a1", "est.a2", "est.a3",
                          "est.bx1", "est.bt1", "est.bxt1", "est.bz",
                          "est.bt2", "est.bt3", "est.bxt2",
                          "est.g11", "est.g12", "est.g13",
                          "est.g21", "est.g22", "est.g23",
                          "est.g31", "est.g32", "est.g33")
  vcov.omtm2 <- mdl.omtm2$vcov
  ses.omtm2 <- c(sqrt(diag(vcov.omtm2)),
                 vcov.omtm2[4, 6], vcov.omtm2[4, 10], vcov.omtm2[6, 10])
  names(ses.omtm2) <- c("se.a1", "se.a2", "se.a3",
                        "se.bx1", "se.bt1", "se.bxt1", "se.bz",
                        "se.bt2", "se.bt3", "se.bxt2",
                        "se.g12", "se.g13",
                        "se.g22", "se.g23",
                        "se.g32", "se.g33",
                        "cov.bx1bxt1", "cov.bx1bxt2", "cov.bxt1bxt2")
  vcov.omtm2 <- vcov.omtm2[c(3, 4:7, 9:10), c(3, 4:7, 9:10)]
  # Calculate probabilities
  pX1Z1 <- expit(-sum(coefs.omtm2[paste0("est.", c("a3", "bx1", "bt1", "bxt1", "bt3", "bxt2", "bz"))]))
  pX1Z0 <- expit(-sum(coefs.omtm2[paste0("est.", c("a3", "bx1", "bt1", "bxt1", "bt3", "bxt2"))]))
  pX0Z1 <- expit(-sum(coefs.omtm2[paste0("est.", c("a3", "bt1", "bt3", "bz"))]))
  pX0Z0 <- expit(-sum(coefs.omtm2[paste0("est.", c("a3", "bt1", "bt3"))]))
  # Calculate the Jacobian matrix
  J0 <- Jt1 <- Jt3 <- ((1 - pX1Z1) * pX1Z1 * pZt0 - (1 - pX0Z1) * pX0Z1 * pZt0) +
    ((1 - pX1Z0) * pX1Z0 * (1 - pZt0) - (1 - pX0Z0) * pX0Z0 * (1 - pZt0))
  Jx <- Jxt1 <- Jxt2 <- (1 - pX1Z1) * pX1Z1 * pZt0 + (1 - pX1Z0) * pX1Z0 * (1 - pZt0)
  Jz <- (1 - pX1Z1) * pX1Z1 * pZt0 - (1 - pX0Z1) * pX0Z1 * pZt0
  J.omtm2 <- matrix(c(J0, Jx, Jt1, Jxt1, Jz, Jt3, Jxt2), nrow = 1)
  RDest.omtm2 <- (pX1Z1 * pZt0 - pX0Z1 * pZt0) +
    (pX1Z0 * (1 - pZt0) - pX0Z0 * (1 - pZt0))
  RDse.omtm2 <- c(sqrt(J.omtm2 %*% vcov.omtm2 %*% t(J.omtm2)))
  cat(paste("OMTM2 convergence code =", code.omtm2, "\n"))

  ##### OMTM2-C #####
  cat("Fitting OMTM2-C \n")
  XMat.ppo2.noX <- XMat.dat.noX[, c(1, 1, 2)]
  params.omtm2.noX <- c(-rev(alpha), -(beta[2] + beta.ppo[2]), -beta[3], -beta[4],
                        beta.ppo[2] - beta.ppo[1], beta.ppo[2], 0,
                        rev(c(t(gamma.mat[[1]]))))
  stepmax <- 1
  code <- 5
  while(code == 5){
    cat(paste("Try stepmax =", stepmax, "\n"))
    reset_step <- FALSE
    tryCatch(
      mdl.omtm2.noX <- omtm1(
        params = params.omtm2.noX, yval = Yomtm2.dat, XMat = XMat.dat.noX, id = id.dat,
        ppo.k = c(2, 3, 3), x.ppo = XMat.ppo2.noX, u = NULL, ref.muc = 2,
        TransIndMtx = NA, ProfileCol = NA, print.level = 1,
        stepmax = stepmax, UseGrad = TRUE, iterlim = 250, check.analyticals = FALSE),
      error = function(e){reset_step <<- TRUE})
    cat(paste("Reset step", reset_step, "\n"))
    if(reset_step == TRUE){
      stepmax <- stepmax * 0.25
    } else {
      stepmax <- stepmax * 1.5
      params.omtm2.noX <- c(mdl.omtm2.noX$alpha, mdl.omtm2.noX$beta, c(t(mdl.omtm2.noX$gamma[[1]])))
      code <- code.omtm2.noX <- mdl.omtm2.noX$control["convergence_code"]
    }}
  coefs.omtm2.noX <- c(mdl.omtm2.noX$alpha, mdl.omtm2.noX$beta, c(t(mdl.omtm2.noX$gamma[[1]])))
  names(coefs.omtm2.noX) <- c("est.a1", "est.a2", "est.a3",
                              "est.bt1", "est.bxt1", "est.bz",
                              "est.bt2", "est.bt3", "est.bxt2",
                              "est.g11", "est.g12", "est.g13",
                              "est.g21", "est.g22", "est.g23",
                              "est.g31", "est.g32", "est.g33")
  vcov.omtm2.noX <- mdl.omtm2.noX$vcov
  ses.omtm2.noX <- c(sqrt(diag(vcov.omtm2.noX)), vcov.omtm2.noX[5, 9])
  names(ses.omtm2.noX) <- c("se.a1", "se.a2", "se.a3",
                            "se.bt1", "se.bxt1", "se.bz",
                            "se.bt2", "se.bt3", "se.bxt2",
                            "se.g12", "se.g13",
                            "se.g22", "se.g23",
                            "se.g32", "se.g33", "cov.bxt1bxt2")
  vcov.omtm2.noX <- vcov.omtm2.noX[c(3, 4:6, 8:9), c(3, 4:6, 8:9)]
  # Calculate probabilities
  pX1Z1 <- expit(-sum(coefs.omtm2.noX[paste0("est.", c("a3", "bt1", "bxt1", "bt3", "bxt2", "bz"))]))
  pX1Z0 <- expit(-sum(coefs.omtm2.noX[paste0("est.", c("a3", "bt1", "bxt1", "bt3", "bxt2"))]))
  pX0Z1 <- expit(-sum(coefs.omtm2.noX[paste0("est.", c("a3", "bt1", "bt3", "bz"))]))
  pX0Z0 <- expit(-sum(coefs.omtm2.noX[paste0("est.", c("a3", "bt1", "bt3"))]))
  # Calculate the Jacobian matrix
  J0 <- Jt1 <- Jt3 <- ((1 - pX1Z1) * pX1Z1 * pZt0 - (1 - pX0Z1) * pX0Z1 * pZt0) +
    ((1 - pX1Z0) * pX1Z0 * (1 - pZt0) - (1 - pX0Z0) * pX0Z0 * (1 - pZt0))
  Jxt1 <- Jxt2 <- (1 - pX1Z1) * pX1Z1 * pZt0 + (1 - pX1Z0) * pX1Z0 * (1 - pZt0)
  Jz <- (1 - pX1Z1) * pX1Z1 * pZt0 - (1 - pX0Z1) * pX0Z1 * pZt0
  J.omtm2.noX <- matrix(c(J0, Jt1, Jxt1, Jz, Jt3, Jxt2), nrow = 1)
  RDest.omtm2.noX <- (pX1Z1 * pZt0 - pX0Z1 * pZt0) +
    (pX1Z0 * (1 - pZt0) - pX0Z0 * (1 - pZt0))
  RDse.omtm2.noX <- c(sqrt(J.omtm2.noX %*% vcov.omtm2.noX %*% t(J.omtm2.noX)))
  cat(paste("OMTM2-C convergence code =", code.omtm2, "\n"))

  ##### OMTM3 #####
  cat("Fitting OMTM3\n")
  XMat.ppo3 <- XMat.dat[, c(2, 2, 3, 3)]
  params.omtm3 <- c(alpha, beta, beta.ppo[1:2], 0, 0, c(t(gamma.mat[[1]])))
  stepmax <- 1
  code <- 5
  while(code == 5){
    cat(paste("Try stepmax =", stepmax, "\n"))
    reset_step <- FALSE
    tryCatch(
      mdl.omtm3 <- omtm1(
        params = params.omtm3, yval = Y.dat, XMat = XMat.dat, id = id.dat,
        ppo.k = c(2, 3, 2, 3), x.ppo = XMat.ppo3, u = NULL, ref.muc = NA,
        TransIndMtx = NA, ProfileCol = NA, print.level = 1,
        stepmax = stepmax, UseGrad = TRUE, iterlim = 250, check.analyticals = FALSE),
      error = function(e){reset_step <<- TRUE})
    cat(paste("Reset step", reset_step, "\n"))
    if(reset_step == TRUE){
      stepmax <- stepmax * 0.25
    } else {
      stepmax <- stepmax * 1.5
      params.omtm3 <- c(mdl.omtm3$alpha, mdl.omtm3$beta, c(t(mdl.omtm3$gamma[[1]])))
      code <- code.omtm3 <- mdl.omtm3$control["convergence_code"]
    }}
  coefs.omtm3 <- c(mdl.omtm3$alpha, mdl.omtm3$beta, c(t(mdl.omtm3$gamma[[1]])))
  names(coefs.omtm3) <- c("est.a1", "est.a2", "est.a3",
                          "est.bx1", "est.bt1", "est.bxt1", "est.bz",
                          "est.bt2", "est.bt3", "est.bxt2", "est.bxt3",
                          "est.g11", "est.g12", "est.g13",
                          "est.g21", "est.g22", "est.g23",
                          "est.g31", "est.g32", "est.g33")
  vcov.omtm3 <- mdl.omtm3$vcov
  ses.omtm3 <- c(sqrt(diag(vcov.omtm3)), vcov.omtm3[4, 6])
  names(ses.omtm3) <- c("se.a1", "se.a2", "se.a3",
                        "se.bx1", "se.bt1", "se.bxt1", "se.bz",
                        "se.bt2", "se.bt3", "se.bxt2", "se.bxt3",
                        "se.g12", "se.g13",
                        "se.g22", "se.g23",
                        "se.g32", "se.g33", "cov.bx1bxt1")
  vcov.omtm3 <- vcov.omtm3[c(1, 4:7), c(1, 4:7)]
  # Calculate probabilities
  pX1Z1 <- expit(sum(coefs.omtm3[paste0("est.", c("a1", "bx1", "bt1", "bxt1", "bz"))]))
  pX1Z0 <- expit(sum(coefs.omtm3[paste0("est.", c("a1", "bx1", "bt1", "bxt1"))]))
  pX0Z1 <- expit(sum(coefs.omtm3[paste0("est.", c("a1", "bt1", "bz"))]))
  pX0Z0 <- expit(sum(coefs.omtm3[paste0("est.", c("a1", "bt1"))]))
  # Calculate the Jacobian matrix
  J0 <- Jt <- ((1 - pX1Z1) * pX1Z1 * pZt0 - (1 - pX0Z1) * pX0Z1 * pZt0) +
    ((1 - pX1Z0) * pX1Z0 * (1 - pZt0) - (1 - pX0Z0) * pX0Z0 * (1 - pZt0))
  Jx <- Jxt <- (1 - pX1Z1) * pX1Z1 * pZt0 + (1 - pX1Z0) * pX1Z0 * (1 - pZt0)
  Jz <- (1 - pX1Z1) * pX1Z1 * pZt0 - (1 - pX0Z1) * pX0Z1 * pZt0
  J.omtm3 <- matrix(c(J0, Jx, Jt, Jxt, Jz), nrow = 1)
  RDest.omtm3 <- (pX1Z1 * pZt0 - pX0Z1 * pZt0) +
    (pX1Z0 * (1 - pZt0) - pX0Z0 * (1 - pZt0))
  RDse.omtm3 <- c(sqrt(J.omtm3 %*% vcov.omtm3 %*% t(J.omtm3)))
  cat(paste("OMTM3 convergence code =", code.omtm3, "\n"))

  ##### OMTM3-C #####
  cat("Fitting OMTM3-C\n")
  XMat.ppo3.noX <- XMat.dat.noX[, c(1, 1, 2, 2)]
  params.omtm3.noX <- c(alpha, beta[-1], beta.ppo[1:2], 0, 0, c(t(gamma.mat[[1]])))
  stepmax <- 1
  code <- 5
  while(code == 5){
    cat(paste("Try stepmax =", stepmax, "\n"))
    reset_step <- FALSE
    tryCatch(
      mdl.omtm3.noX <- omtm1(
        params = params.omtm3.noX, yval = Y.dat, XMat = XMat.dat.noX, id = id.dat,
        ppo.k = c(2, 3, 2, 3), x.ppo = XMat.ppo3.noX, u = NULL, ref.muc = NA,
        TransIndMtx = NA, ProfileCol = NA, print.level = 1,
        stepmax = stepmax, UseGrad = TRUE, iterlim = 250, check.analyticals = FALSE),
      error = function(e){reset_step <<- TRUE})
    cat(paste("Reset step", reset_step, "\n"))
    if(reset_step == TRUE){
      stepmax <- stepmax * 0.25
    } else {
      stepmax <- stepmax * 1.5
      params.omtm3.noX <- c(mdl.omtm3.noX$alpha, mdl.omtm3.noX$beta, c(t(mdl.omtm3.noX$gamma[[1]])))
      code <- code.omtm3.noX <- mdl.omtm3.noX$control["convergence_code"]
    }}
  coefs.omtm3.noX <- c(mdl.omtm3.noX$alpha, mdl.omtm3.noX$beta, c(t(mdl.omtm3.noX$gamma[[1]])))
  names(coefs.omtm3.noX) <- c("est.a1", "est.a2", "est.a3",
                              "est.bt1", "est.bxt1", "est.bz",
                              "est.bt2", "est.bt3", "est.bxt2", "est.bxt3",
                              "est.g11", "est.g12", "est.g13",
                              "est.g21", "est.g22", "est.g23",
                              "est.g31", "est.g32", "est.g33")
  vcov.omtm3.noX <- mdl.omtm3.noX$vcov
  ses.omtm3.noX <- sqrt(diag(vcov.omtm3.noX))
  names(ses.omtm3.noX) <- c("se.a1", "se.a2", "se.a3",
                            "se.bt1", "se.bxt1", "se.bz",
                            "se.bt2", "se.bt3", "se.bxt2", "se.bxt3",
                            "se.g12", "se.g13",
                            "se.g22", "se.g23",
                            "se.g32", "se.g33")
  vcov.omtm3.noX <- vcov.omtm3.noX[c(1, 4:6), c(1, 4:6)]
  # Calculate probabilities
  pX1Z1 <- expit(sum(coefs.omtm3.noX[paste0("est.", c("a1", "bt1", "bxt1", "bz"))]))
  pX1Z0 <- expit(sum(coefs.omtm3.noX[paste0("est.", c("a1", "bt1", "bxt1"))]))
  pX0Z1 <- expit(sum(coefs.omtm3.noX[paste0("est.", c("a1", "bt1", "bz"))]))
  pX0Z0 <- expit(sum(coefs.omtm3.noX[paste0("est.", c("a1", "bt1"))]))
  # Calculate the Jacobian matrix
  J0 <- Jt <- ((1 - pX1Z1) * pX1Z1 * pZt0 - (1 - pX0Z1) * pX0Z1 * pZt0) +
    ((1 - pX1Z0) * pX1Z0 * (1 - pZt0) - (1 - pX0Z0) * pX0Z0 * (1 - pZt0))
  Jxt <- (1 - pX1Z1) * pX1Z1 * pZt0 + (1 - pX1Z0) * pX1Z0 * (1 - pZt0)
  Jz <- (1 - pX1Z1) * pX1Z1 * pZt0 - (1 - pX0Z1) * pX0Z1 * pZt0
  J.omtm3.noX <- matrix(c(J0, Jt, Jxt, Jz), nrow = 1)
  RDest.omtm3.noX <- (pX1Z1 * pZt0 - pX0Z1 * pZt0) +
    (pX1Z0 * (1 - pZt0) - pX0Z0 * (1 - pZt0))
  RDse.omtm3.noX <- c(sqrt(J.omtm3.noX %*% vcov.omtm3.noX %*% t(J.omtm3.noX)))
  cat(paste("OMTM3-C convergence code =", code.omtm1.noX, "\n"))

  ##### Survival ####
  cat("Fitting Cox PH\n")
  # Fit Cox PH model and calculate exp(eta)
  mdl.cox <- coxph(Surv(time, Ybin) ~ X + Z, data = dat.cox, ties = "breslow")
  new.dat <- expand.grid(X = c(0, 1), Z = c(0, 1))
  surv.breslow <- survfit(mdl.cox, newdata = new.dat)
  beta_hat <- coef(mdl.cox)
  dat.cox$exp_eta <- exp(c(t(beta_hat) %*% t(dat.cox[,c("X", "Z")])))
  lambda0_hat <- sapply(0:27, function (t) {
    dat.tmp1 <- dat.cox %>% filter(time == t & Ybin == 1)
    dat.tmp0 <- dat.cox %>% filter(time >= t)
    nrow(dat.tmp1) / sum(dat.tmp0$exp_eta)
  })
  # Breslow estimator for cumulative baseline hazard
  Lambda0_hat <- cumsum(lambda0_hat)
  # Variance
  # for beta_hat
  beta_hat.var <- vcov(mdl.cox)
  # for Lambda0_hat
  d <- c(); tmp0 <- c(); tmp1x <- c(); tmp1z <- c()
  for (t in 0:27){
    dat.tmp1 <- dat.cox %>% filter(time == t & Ybin == 1)
    dat.tmp0 <- dat.cox %>% filter(time >= t)
    tmp0 <- c(tmp0, sum(dat.tmp0$exp_eta))
    tmp1x <- c(tmp1x, sum(dat.tmp0$X * dat.tmp0$exp_eta))
    tmp1z <- c(tmp1z, sum(dat.tmp0$Z * dat.tmp0$exp_eta))
    d <- c(d, nrow(dat.tmp1))
  }
  Lambda0_hat.var <- c(); cov_hat <- matrix(nrow=mi, ncol=2)
  for (t in 0:27) {
    mat.tmp <- matrix(c(sum(tmp1x[1:(t+1)] / tmp0[1:(t+1)]^2 * d[1:(t+1)]),
                        sum(tmp1z[1:(t+1)] / tmp0[1:(t+1)]^2 * d[1:(t+1)])), nrow = 1)
    Lambda0_hat.var <- c(Lambda0_hat.var,
                         c(mat.tmp %*% beta_hat.var %*% t(mat.tmp)) + sum(d[1:(t+1)] / tmp0[1:(t+1)]^2))
    # for Cov(Lambda0_hat, beta_hat)
    cov_hat[t+1, ] <- - mat.tmp %*% beta_hat.var
  }
  # Calculate cumulative hazard function at t
  t <- 27
  LambdaX0Z0 <- Lambda0_hat[t+1]
  LambdaX1Z1 <- LambdaX0Z0 * exp(sum(beta_hat))
  LambdaX1Z0 <- LambdaX0Z0 * exp(beta_hat[1])
  LambdaX0Z1 <- LambdaX0Z0 * exp(beta_hat[2])
  # Calculate survival function at t
  SX0Z0 <- exp(-LambdaX0Z0)
  SX1Z1 <- exp(-LambdaX1Z1)
  SX1Z0 <- exp(-LambdaX1Z0)
  SX0Z1 <- exp(-LambdaX0Z1)
  # RD estimation
  RDest.surv <- (SX0Z1 * pZt0 - SX1Z1 * pZt0) + (SX0Z0 * (1 - pZt0) - SX1Z0 * (1 - pZt0))
  # Covariance matrix
  Sigmat_hat <- rbind(matrix(c(Lambda0_hat.var[t+1], cov_hat[t+1, ]), nrow = 1),
                      cbind(t(cov_hat[t+1, , drop = F]), beta_hat.var))
  # Jacobian matrix
  Jt <- matrix(c(
    (-exp(beta_hat[2]) * SX0Z1 * pZt0 + exp(sum(beta_hat)) * SX1Z1 * pZt0) +
      (-SX0Z0 * (1 - pZt0) + exp(beta_hat[1]) * SX1Z0 * (1 - pZt0)),
    LambdaX0Z0 * (exp(sum(beta_hat)) * SX1Z1 * pZt0 + exp(beta_hat[1]) * SX1Z0 * (1 - pZt0)),
    LambdaX0Z0 * (exp(sum(beta_hat)) * SX1Z1 * pZt0 - exp(beta_hat[2]) * SX0Z1 * pZt0)), nrow = 1)
  RDse.surv <- c(sqrt(Jt %*% Sigmat_hat %*% t(Jt)))

  ##### Save #####
  cat("Saving Results\n")
  # log(OR)
  coefs.rst <- bind_rows(
    coefs.glmb, coefs.bmtm, coefs.bmtm.noX, coefs.glmo1, coefs.glmo2, coefs.glmo3,
    coefs.omtm1, coefs.omtm1.noX, coefs.omtm2, coefs.omtm2.noX, coefs.omtm3, coefs.omtm3.noX)
  ses.rst <- bind_rows(
    ses.glmb, ses.bmtm, ses.bmtm.noX, ses.glmo1, ses.glmo2, ses.glmo3,
    ses.omtm1, ses.omtm1.noX, ses.omtm2, ses.omtm2.noX, ses.omtm3, ses.omtm3.noX)
  sim1rst <- data.frame(
    data_gen = "PPO T and XT",
    fu       = "Complete",
    mdl_fit  = c("GLMB", "BMTM", "BMTM-C", "GLMO1", "GLMO2", "GLMO3",
                 "OMTM1", "OMTM1-C", "OMTM2", "OMTM2-C", "OMTM3", "OMTM3-C")
  )
  sim1rst <- cbind(sim1rst, coefs.rst, ses.rst)
  if (file.exists("output-missX")){
    save(sim1rst, file = paste0("output-missX/sim1-run-", run, ".RData"))}
  # RD
  sim2rst <- data.frame(
    data_gen = "PPO T and XT",
    fu       = "Complete",
    mdl_fit  = c("GLMB", "BMTM", "BMTM-C", "GLMO1", "GLMO2", "GLMO3",
                 "OMTM1", "OMTM1-C", "OMTM2", "OMTM2-C", "OMTM3", "OMTM3-C", "Survival"),
    RDest    = c(RDest.glmb, RDest.bmtm, RDest.bmtm.noX, RDest.glmo1, RDest.glmo2, RDest.glmo3,
                 RDest.omtm1, RDest.omtm1.noX, RDest.omtm2, RDest.omtm2.noX, RDest.omtm3, RDest.omtm3.noX, RDest.surv),
    RDse     = c(RDse.glmb, RDse.bmtm, RDse.bmtm.noX, RDse.glmo1, RDse.glmo2, RDse.glmo3,
                 RDse.omtm1, RDse.omtm1.noX, RDse.omtm2, RDse.omtm2.noX, RDse.omtm3, RDse.omtm3.noX, RDse.surv),
    code     = c(NA, code.bmtm, code.bmtm.noX, rep(NA, 3),
                 code.omtm1, code.omtm1.noX, code.omtm2, code.omtm2.noX, code.omtm3, code.omtm3.noX, NA),
    pZt0     = pZt0
  )
  if (file.exists("output-missX")){
    save(sim2rst, file = paste0("output-missX/sim2-run-", run, ".RData"))}
  
  ##### Data - Complete (misspecification of beta_z) ####
  cat("Generating random data from known parameters\n")
  # Generate binary treatment indicator X (`tx`)
  tx <- rep(rbinom(N, size = 1, prob = pX), each = mi)
  # Generate binary baseline status Z (Fix p(Z) = 0.5)
  Z <- rep(0, N)
  Z[sample(1:N, N * pZ)] <- 1
  Z <- rep(Z, each = mi)
  # Construct the design matrix `XMat` in the form of (tx, logt, tx*logt, Z)
  XMat <- cbind(tx, logt, tx*logt, Z)
  # Generate outcome `Y` by relaxing the PO assumption on T and Z
  # in the presence of the absorbing state (k=1)
  XMat.ppo <- XMat[, c(2, 2, 4, 4)]
  beta.ppo <- c(-1.75, -4, 0, -1)
  ppo.k <- c(2, 3, 2, 3)
  id <- rep(1:N, each = mi)
  Y <- GenDatOMTM1.ppo(
    id = id, XMat = XMat, alpha = alpha, beta = beta, gamma.mat.list = gamma.mat,
    ppo.k = ppo.k, XMat.ppo = XMat.ppo, beta.ppo = beta.ppo)
  # Incorporate the generated data into a dataframe
  dat <- data.frame(id = id, X = tx, time = t, logt = logt, Z = Z, Y = Y)
  # Create binary outcome variable: Ybin=1 (event Y=1), Ybin=0 (no event Y=2:4)
  dat$Ybin <- ifelse(dat$Y == 1, 1, 0)
  dat0 <- dat
  id.event <- dat[with(dat, logt == 1 & Ybin == 1), ]$id
  dat.noevent <- dat %>% filter(!(id %in% id.event)) %>%
    group_by(id) %>% filter(time == max(time)) %>% ungroup()
  dat.event <- dat %>% filter(id %in% id.event & Ybin == 1) %>%
    group_by(id) %>% filter(time == min(time)) %>% ungroup()
  dat.cox <- rbind(dat.noevent, dat.event)
  # Save data
  if (file.exists("data-missZ")){
    save(dat, file = paste0("data-missZ/data-", run, ".RData"))}
  
  # Calculate p(Z) at the beginning of the study
  dat.t0 <- dat %>% filter(logt == 0)
  pZt0 <- mean(dat.t0$Z == 1)
  
  ##### GLMB ####
  cat("Fitting GLMB\n")
  # Subset data at the end of follow-up time only
  dat.cross <- subset(dat, logt == 1)
  dat.cross$Ybin <- ifelse(dat.cross$Ybin == 1, 1, 0)
  # Fit GLMB model and extract coefficients, covariance matrix
  mdl.glmb <- glm(Ybin ~ X + Z, data = dat.cross, family = "binomial")
  coefs.glmb <- coef(mdl.glmb)  # (beta_0, beta_x, beta_z)
  names(coefs.glmb) <- c("est.a1", "est.bx1", "est.bz")
  vcov.glmb <- vcov(mdl.glmb)
  ses.glmb <- sqrt(diag(vcov.glmb))
  names(ses.glmb) <- c("se.a1", "se.bx1", "se.bz")
  # Calculate probabilities
  pX1Z1 <- expit(sum(coefs.glmb))
  pX1Z0 <- expit(sum(coefs.glmb[c(1, 2)]))
  pX0Z1 <- expit(sum(coefs.glmb[c(1, 3)]))
  pX0Z0 <- expit(sum(coefs.glmb[1]))
  # Calculate the Jacobian matrix
  J.glmb <- matrix(c(
    ((1 - pX1Z1) * pX1Z1 * pZt0 - (1 - pX0Z1) * pX0Z1 * pZt0) +
      ((1 - pX1Z0) * pX1Z0 * (1 - pZt0) - (1 - pX0Z0) * pX0Z0 * (1 - pZt0)),
    (1 - pX1Z1) * pX1Z1 * pZt0 + (1 - pX1Z0) * pX1Z0 * (1 - pZt0),
    (1 - pX1Z1) * pX1Z1 * pZt0 - (1 - pX0Z1) * pX0Z1 * pZt0), nrow = 1)
  # Calculate RD estimate and its standard error
  RDest.glmb <- (pX1Z1 * pZt0 - pX0Z1 * pZt0) +
    (pX1Z0 * (1 - pZt0) - pX0Z0 * (1 - pZt0))
  RDse.glmb <- c(sqrt(J.glmb %*% vcov.glmb %*% t(J.glmb)))
  
  ##### BMTM ####
  cat("Fitting BMTM\n")
  # Create binary outcome variable: Ybin=1 (event Y=1), Ybin=2 (no event Y=2:4)
  dat$Ybin <- ifelse(dat$Y != 1, 2, 1)
  # Fit BMTM model and extract coefficients, covariance matrix
  XMat.dat <- with(dat, cbind(X, logt, X*logt, Z))
  Ybin.dat <- dat$Ybin
  id.dat <- dat$id
  params.bmtm <- c(alpha[1], beta, gamma.mat[[1]][1, 1])
  stepmax <- 1
  code <- 5
  while(code == 5){
    cat(paste("Try stepmax =", stepmax, "\n"))
    reset_step <- FALSE
    tryCatch(
      mdl.bmtm <- omtm1(
        params = params.bmtm, yval = Ybin.dat, XMat = XMat.dat, id = id.dat,
        ppo.k = NULL, x.ppo = NULL, u = NULL, ref.muc = NA,
        TransIndMtx = NA, ProfileCol = NA, print.level = 1,
        stepmax = stepmax, UseGrad = TRUE, iterlim = 250, check.analyticals = FALSE),
      error = function(e){reset_step <<- TRUE})
    cat(paste("Reset step", reset_step, "\n"))
    if(reset_step == TRUE){
      stepmax <- stepmax * 0.25
    } else {
      stepmax <- stepmax * 1.5
      params.bmtm <- c(mdl.bmtm$alpha, mdl.bmtm$beta, mdl.bmtm$gamma[[1]][1, 1])
      code <- code.bmtm <- mdl.bmtm$control["convergence_code"]
    }}
  coefs.bmtm <- c(mdl.bmtm$alpha, mdl.bmtm$beta, mdl.bmtm$gamma[[1]][1, 1])
  names(coefs.bmtm) <- c("est.a1", "est.bx1", "est.bt1", "est.bxt1", "est.bz", "est.g11")
  vcov.bmtm <- mdl.bmtm$vcov
  ses.bmtm <- c(sqrt(diag(vcov.bmtm)), vcov.bmtm[2, 4])
  names(ses.bmtm) <- c("se.a1", "se.bx1", "se.bt1", "se.bxt1", "se.bz", "cov.bx1bxt1")
  # Calculate probabilities
  pX1Z1 <- expit(sum(coefs.bmtm[paste0("est.", c("a1", "bx1", "bt1", "bxt1", "bz"))]))
  pX1Z0 <- expit(sum(coefs.bmtm[paste0("est.", c("a1", "bx1", "bt1", "bxt1"))]))
  pX0Z1 <- expit(sum(coefs.bmtm[paste0("est.", c("a1", "bt1", "bz"))]))
  pX0Z0 <- expit(sum(coefs.bmtm[paste0("est.", c("a1", "bt1"))]))
  # Calculate the Jacobian matrix
  J0 <- Jt <- ((1 - pX1Z1) * pX1Z1 * pZt0 - (1 - pX0Z1) * pX0Z1 * pZt0) +
    ((1 - pX1Z0) * pX1Z0 * (1 - pZt0) - (1 - pX0Z0) * pX0Z0 * (1 - pZt0))
  Jx <- Jxt <- (1 - pX1Z1) * pX1Z1 * pZt0 + (1 - pX1Z0) * pX1Z0 * (1 - pZt0)
  Jz <- (1 - pX1Z1) * pX1Z1 * pZt0 - (1 - pX0Z1) * pX0Z1 * pZt0
  J.bmtm <- matrix(c(J0, Jx, Jt, Jxt, Jz), nrow = 1)
  RDest.bmtm <- (pX1Z1 * pZt0 - pX0Z1 * pZt0) +
    (pX1Z0 * (1 - pZt0) - pX0Z0 * (1 - pZt0))
  RDse.bmtm <- c(sqrt(J.bmtm %*% vcov.bmtm %*% t(J.bmtm)))
  cat(paste("BMTM convergence code =", code.bmtm, "\n"))
  
  ##### BMTM-C ####
  cat("Fitting BMTM-C\n")
  # Fit BMTM-C model and extract coefficients, covariance matrix
  XMat.dat.noX <- with(dat, cbind(logt, X*logt, Z))
  params.bmtm.noX <- c(alpha[1], beta[-1], gamma.mat[[1]][1, 1])
  stepmax <- 1
  code <- 5
  while(code == 5){
    cat(paste("Try stepmax =", stepmax, "\n"))
    reset_step <- FALSE
    tryCatch(
      mdl.bmtm.noX <- omtm1(
        params = params.bmtm.noX, yval = Ybin.dat, XMat = XMat.dat.noX, id = id.dat,
        ppo.k = NULL, x.ppo = NULL, u = NULL, ref.muc = NA,
        TransIndMtx = NA, ProfileCol = NA, print.level = 1,
        stepmax = 1, UseGrad = TRUE, iterlim = 250, check.analyticals = FALSE),
      error = function(e){reset_step <<- TRUE})
    cat(paste("Reset step", reset_step, "\n"))
    if(reset_step == TRUE){
      stepmax <- stepmax * 0.25
    } else {
      stepmax <- stepmax * 1.5
      params.bmtm.noX <- c(mdl.bmtm.noX$alpha, mdl.bmtm.noX$beta, mdl.bmtm.noX$gamma[[1]][1, 1])
      code <- code.bmtm.noX <- mdl.bmtm.noX$control["convergence_code"]
    }}
  coefs.bmtm.noX <- c(mdl.bmtm.noX$alpha, mdl.bmtm.noX$beta, mdl.bmtm.noX$gamma[[1]][1, 1])
  names(coefs.bmtm.noX) <- c("est.a1", "est.bt1", "est.bxt1", "est.bz", "est.g11")
  vcov.bmtm.noX <- mdl.bmtm.noX$vcov
  ses.bmtm.noX <- sqrt(diag(vcov.bmtm.noX))
  names(ses.bmtm.noX) <- c("se.a1", "se.bt1", "se.bxt1", "se.bz")
  # Calculate probabilities
  pX1Z1 <- expit(sum(coefs.bmtm.noX[paste0("est.", c("a1", "bt1", "bxt1", "bz"))]))
  pX1Z0 <- expit(sum(coefs.bmtm.noX[paste0("est.", c("a1", "bt1", "bxt1"))]))
  pX0Z1 <- expit(sum(coefs.bmtm.noX[paste0("est.", c("a1", "bt1", "bz"))]))
  pX0Z0 <- expit(sum(coefs.bmtm.noX[paste0("est.", c("a1", "bt1"))]))
  # Calculate the Jacobian matrix
  J0 <- Jt <- ((1 - pX1Z1) * pX1Z1 * pZt0 - (1 - pX0Z1) * pX0Z1 * pZt0) +
    ((1 - pX1Z0) * pX1Z0 * (1 - pZt0) - (1 - pX0Z0) * pX0Z0 * (1 - pZt0))
  Jxt <- (1 - pX1Z1) * pX1Z1 * pZt0 + (1 - pX1Z0) * pX1Z0 * (1 - pZt0)
  Jz <- (1 - pX1Z1) * pX1Z1 * pZt0 - (1 - pX0Z1) * pX0Z1 * pZt0
  J.bmtm.noX <- matrix(c(J0, Jt, Jxt, Jz), nrow=1)
  RDest.bmtm.noX <- (pX1Z1 * pZt0 - pX0Z1 * pZt0) +
    (pX1Z0 * (1 - pZt0) - pX0Z0 * (1 - pZt0))
  RDse.bmtm.noX <- c(sqrt(J.bmtm.noX %*% vcov.bmtm.noX %*% t(J.bmtm.noX)))
  cat(paste("BMTM-C convergence code =", code.bmtm.noX, "\n"))
  
  ##### GLMO1 #####
  # Fit GLMB model and extract coefficients, covariance matrix
  mdl.glmo1 <- vglm(ordered(Y, levels = 1:4) ~ X + Z, data = dat.cross,
                    family = cumulative(parallel = TRUE))
  coefs.glmo1 <- coef(mdl.glmo1)
  names(coefs.glmo1) <- c("est.a1", "est.a2", "est.a3", "est.bx1", "est.bz")
  vcov.glmo1 <- vcov(mdl.glmo1)
  ses.glmo1 <- sqrt(diag(vcov.glmo1))
  names(ses.glmo1) <- c("se.a1", "se.a2", "se.a3", "se.bx1", "se.bz")
  vcov.glmo1 <- vcov.glmo1[c(1, 4:5),c(1, 4:5)]
  # Calculate probabilities
  pX1Z1 <- expit(sum(coefs.glmo1[paste0("est.", c("a1", "bx1", "bz"))]))
  pX1Z0 <- expit(sum(coefs.glmo1[paste0("est.", c("a1", "bx1"))]))
  pX0Z1 <- expit(sum(coefs.glmo1[paste0("est.", c("a1", "bz"))]))
  pX0Z0 <- expit(sum(coefs.glmo1[paste0("est.", c("a1"))]))
  # Calculate the Jacobian matrix
  J.glmo1 <- matrix(c(
    ((1 - pX1Z1) * pX1Z1 * pZt0 - (1 - pX0Z1) * pX0Z1 * pZt0) +
      ((1 - pX1Z0) * pX1Z0 * (1 - pZt0) - (1 - pX0Z0) * pX0Z0 * (1 - pZt0)),
    (1 - pX1Z1) * pX1Z1 * pZt0 + (1 - pX1Z0) * pX1Z0 * (1 - pZt0),
    (1 - pX1Z1) * pX1Z1 * pZt0 - (1 - pX0Z1) * pX0Z1 * pZt0), nrow = 1)
  # Calculate RD estimate and its standard error
  RDest.glmo1 <- (pX1Z1 * pZt0 - pX0Z1 * pZt0) +
    (pX1Z0 * (1 - pZt0) - pX0Z0 * (1 - pZt0))
  RDse.glmo1 <- c(sqrt(J.glmo1 %*% vcov.glmo1 %*% t(J.glmo1)))
  
  ##### GLMO2 #####
  # Fit GLMB model and extract coefficients, covariance matrix
  mat <- cbind(rbind(1, 0, 0), rbind(0, 1, 1))
  clist <- list("(Intercept)" = diag(3), "X" = mat, "Z" = rbind(1, 1, 1))
  mdl.glmo2 <- vglm(ordered(Y, levels = 1:4) ~ X + Z, data = dat.cross,
                    constraints = clist, family = cumulative(parallel = FALSE))
  coefs.glmo2 <- coef(mdl.glmo2)
  names(coefs.glmo2) <- c("est.a1", "est.a2", "est.a3", "est.bx1", "est.bx2", "est.bz")
  vcov.glmo2 <- vcov(mdl.glmo2)
  ses.glmo2 <- sqrt(diag(vcov.glmo2))
  names(ses.glmo2) <- c("se.a1", "se.a2", "se.a3", "se.bx1", "se.bx2", "se.bz")
  vcov.glmo2 <- vcov.glmo2[c(1, 4, 6), c(1, 4, 6)]
  # Calculate probabilities
  pX1Z1 <- expit(sum(coefs.glmo2[paste0("est.", c("a1", "bx1", "bz"))]))
  pX1Z0 <- expit(sum(coefs.glmo2[paste0("est.", c("a1", "bx1"))]))
  pX0Z1 <- expit(sum(coefs.glmo2[paste0("est.", c("a1", "bz"))]))
  pX0Z0 <- expit(sum(coefs.glmo2[paste0("est.", c("a1"))]))
  # Calculate the Jacobian matrix
  J.glmo2 <- matrix(c(
    ((1 - pX1Z1) * pX1Z1 * pZt0 - (1 - pX0Z1) * pX0Z1 * pZt0) +
      ((1 - pX1Z0) * pX1Z0 * (1 - pZt0) - (1 - pX0Z0) * pX0Z0 * (1 - pZt0)),
    (1 - pX1Z1) * pX1Z1 * pZt0 + (1 - pX1Z0) * pX1Z0 * (1 - pZt0),
    (1 - pX1Z1) * pX1Z1 * pZt0 - (1 - pX0Z1) * pX0Z1 * pZt0), nrow = 1)
  # Calculate RD estimate and its standard error
  RDest.glmo2 <- (pX1Z1 * pZt0 - pX0Z1 * pZt0) +
    (pX1Z0 * (1 - pZt0) - pX0Z0 * (1 - pZt0))
  RDse.glmo2 <- c(sqrt(J.glmo2 %*% vcov.glmo2 %*% t(J.glmo2)))
  
  ##### GLMO3 #####
  # Fit GLMB model and extract coefficients, covariance matrix
  mdl.glmo3 <- vglm(ordered(Y, levels = 1:4) ~ X + Z, data = dat.cross,
                    family = cumulative(parallel = FALSE ~ X))
  coefs.glmo3 <- coef(mdl.glmo3)
  names(coefs.glmo3) <- c("est.a1", "est.a2", "est.a3", "est.bx1", "est.bx2", "est.bx3", "est.bz")
  vcov.glmo3 <- vcov(mdl.glmo3)
  ses.glmo3 <- sqrt(diag(vcov.glmo3))
  names(ses.glmo3) <- c("se.a1", "se.a2", "se.a3", "se.bx1", "se.bx2", "se.bx3", "se.bz")
  vcov.glmo3 <- vcov.glmo3[c(1, 4, 7), c(1, 4, 7)]
  # Calculate probabilities
  pX1Z1 <- expit(sum(coefs.glmo3[paste0("est.", c("a1", "bx1", "bz"))]))
  pX1Z0 <- expit(sum(coefs.glmo3[paste0("est.", c("a1", "bx1"))]))
  pX0Z1 <- expit(sum(coefs.glmo3[paste0("est.", c("a1", "bz"))]))
  pX0Z0 <- expit(sum(coefs.glmo3[paste0("est.", c("a1"))]))
  # Calculate the Jacobian matrix
  J.glmo3 <- matrix(c(
    ((1 - pX1Z1) * pX1Z1 * pZt0 - (1 - pX0Z1) * pX0Z1 * pZt0) +
      ((1 - pX1Z0) * pX1Z0 * (1 - pZt0) - (1 - pX0Z0) * pX0Z0 * (1 - pZt0)),
    (1 - pX1Z1) * pX1Z1 * pZt0 + (1 - pX1Z0) * pX1Z0 * (1 - pZt0),
    (1 - pX1Z1) * pX1Z1 * pZt0 - (1 - pX0Z1) * pX0Z1 * pZt0), nrow = 1)
  # Calculate RD estimate and its standard error
  RDest.glmo3 <- (pX1Z1 * pZt0 - pX0Z1 * pZt0) +
    (pX1Z0 * (1 - pZt0) - pX0Z0 * (1 - pZt0))
  RDse.glmo3 <- c(sqrt(J.glmo3 %*% vcov.glmo3 %*% t(J.glmo3)))
  
  ##### OMTM1 #####
  cat("Fitting OMTM1\n")
  XMat.ppo1 <- XMat.dat[, c(2, 2)]
  Y.dat <- dat$Y
  params.omtm1 <- c(alpha, beta, beta.ppo[1:2], c(t(gamma.mat[[1]])))
  stepmax <- 1
  code <- 5
  while(code == 5){
    cat(paste("Try stepmax =", stepmax, "\n"))
    reset_step <- FALSE
    tryCatch(
      mdl.omtm1 <- omtm1(
        params = params.omtm1, yval = Y.dat, XMat = XMat.dat, id = id.dat,
        ppo.k = c(2, 3), x.ppo = XMat.ppo1, u = NULL, ref.muc = NA,
        TransIndMtx = NA, ProfileCol = NA, print.level = 1,
        stepmax = stepmax, UseGrad = TRUE, iterlim = 250, check.analyticals = FALSE),
      error = function(e){reset_step <<- TRUE})
    cat(paste("Reset step", reset_step, "\n"))
    if(reset_step == TRUE){
      stepmax <- stepmax * 0.25
    } else {
      stepmax <- stepmax * 1.5
      params.omtm1 <- c(mdl.omtm1$alpha, mdl.omtm1$beta, c(t(mdl.omtm1$gamma[[1]])))
      code <- code.omtm1 <- mdl.omtm1$control["convergence_code"]
    }}
  coefs.omtm1 <- c(mdl.omtm1$alpha, mdl.omtm1$beta, c(t(mdl.omtm1$gamma[[1]])))
  names(coefs.omtm1) <- c("est.a1", "est.a2", "est.a3",
                          "est.bx1", "est.bt1", "est.bxt1", "est.bz",
                          "est.bt2", "est.bt3",
                          "est.g11", "est.g12", "est.g13",
                          "est.g21", "est.g22", "est.g23",
                          "est.g31", "est.g32", "est.g33")
  vcov.omtm1 <- mdl.omtm1$vcov
  ses.omtm1 <- c(sqrt(diag(vcov.omtm1)), vcov.omtm1[4, 6])
  names(ses.omtm1) <- c("se.a1", "se.a2", "se.a3",
                        "se.bx1", "se.bt1", "se.bxt1", "se.bz",
                        "se.bt2", "se.bt3",
                        "se.g12", "se.g13",
                        "se.g22", "se.g23",
                        "se.g32", "se.g33", "cov.bx1bxt1")
  vcov.omtm1 <- vcov.omtm1[c(1, 4:7), c(1, 4:7)]
  # Calculate probabilities
  pX1Z1 <- expit(sum(coefs.omtm1[paste0("est.", c("a1", "bx1", "bt1", "bxt1", "bz"))]))
  pX1Z0 <- expit(sum(coefs.omtm1[paste0("est.", c("a1", "bx1", "bt1", "bxt1"))]))
  pX0Z1 <- expit(sum(coefs.omtm1[paste0("est.", c("a1", "bt1", "bz"))]))
  pX0Z0 <- expit(sum(coefs.omtm1[paste0("est.", c("a1", "bt1"))]))
  # Calculate the Jacobian matrix
  J0 <- Jt <- ((1 - pX1Z1) * pX1Z1 * pZt0 - (1 - pX0Z1) * pX0Z1 * pZt0) +
    ((1 - pX1Z0) * pX1Z0 * (1 - pZt0) - (1 - pX0Z0) * pX0Z0 * (1 - pZt0))
  Jx <- Jxt <- (1 - pX1Z1) * pX1Z1 * pZt0 + (1 - pX1Z0) * pX1Z0 * (1 - pZt0)
  Jz <- (1 - pX1Z1) * pX1Z1 * pZt0 - (1 - pX0Z1) * pX0Z1 * pZt0
  J.omtm1 <- matrix(c(J0, Jx, Jt, Jxt, Jz), nrow = 1)
  RDest.omtm1 <- (pX1Z1 * pZt0 - pX0Z1 * pZt0) +
    (pX1Z0 * (1 - pZt0) - pX0Z0 * (1 - pZt0))
  RDse.omtm1 <- c(sqrt(J.omtm1 %*% vcov.omtm1 %*% t(J.omtm1)))
  cat(paste("OMTM1 convergence code =", code.omtm1, "\n"))
  
  ##### OMTM1-C #####
  cat("Fitting OMTM1-C\n")
  XMat.ppo1.noX <- XMat.dat.noX[, c(1, 1)]
  params.omtm1.noX <- c(alpha, beta[-1], beta.ppo[1:2], c(t(gamma.mat[[1]])))
  stepmax <- 1
  code <- 5
  while(code == 5){
    cat(paste("Try stepmax =", stepmax, "\n"))
    reset_step <- FALSE
    tryCatch(
      mdl.omtm1.noX <- omtm1(
        params = params.omtm1.noX, yval = Y.dat, XMat = XMat.dat.noX, id = id.dat,
        ppo.k = c(2, 3), x.ppo = XMat.ppo1.noX, u = NULL, ref.muc = NA,
        TransIndMtx = NA, ProfileCol = NA, print.level = 1,
        stepmax = stepmax, UseGrad = TRUE, iterlim = 250, check.analyticals = FALSE),
      error = function(e){reset_step <<- TRUE})
    cat(paste("Reset step", reset_step, "\n"))
    if(reset_step == TRUE){
      stepmax <- stepmax * 0.25
    } else {
      stepmax <- stepmax * 1.5
      params.omtm1.noX <- c(mdl.omtm1.noX$alpha, mdl.omtm1.noX$beta, c(t(mdl.omtm1.noX$gamma[[1]])))
      code <- code.omtm1.noX <- mdl.omtm1.noX$control["convergence_code"]
    }}
  coefs.omtm1.noX <- c(mdl.omtm1.noX$alpha, mdl.omtm1.noX$beta, c(t(mdl.omtm1.noX$gamma[[1]])))
  names(coefs.omtm1.noX) <- c("est.a1", "est.a2", "est.a3",
                              "est.bt1", "est.bxt1", "est.bz",
                              "est.bt2", "est.bt3",
                              "est.g11", "est.g12", "est.g13",
                              "est.g21", "est.g22", "est.g23",
                              "est.g31", "est.g32", "est.g33")
  vcov.omtm1.noX <- mdl.omtm1.noX$vcov
  ses.omtm1.noX <- sqrt(diag(vcov.omtm1.noX))
  names(ses.omtm1.noX) <- c("se.a1", "se.a2", "se.a3",
                            "se.bt1", "se.bxt1", "se.bz",
                            "se.bt2", "se.bt3",
                            "se.g12", "se.g13",
                            "se.g22", "se.g23",
                            "se.g32", "se.g33")
  vcov.omtm1.noX <- vcov.omtm1.noX[c(1, 4:6), c(1, 4:6)]
  # Calculate probabilities
  pX1Z1 <- expit(sum(coefs.omtm1.noX[paste0("est.", c("a1", "bt1", "bxt1", "bz"))]))
  pX1Z0 <- expit(sum(coefs.omtm1.noX[paste0("est.", c("a1", "bt1", "bxt1"))]))
  pX0Z1 <- expit(sum(coefs.omtm1.noX[paste0("est.", c("a1", "bt1", "bz"))]))
  pX0Z0 <- expit(sum(coefs.omtm1.noX[paste0("est.", c("a1", "bt1"))]))
  # Calculate the Jacobian matrix
  J0 <- Jt <- ((1 - pX1Z1) * pX1Z1 * pZt0 - (1 - pX0Z1) * pX0Z1 * pZt0) +
    ((1 - pX1Z0) * pX1Z0 * (1 - pZt0) - (1 - pX0Z0) * pX0Z0 * (1 - pZt0))
  Jxt <- (1 - pX1Z1) * pX1Z1 * pZt0 + (1 - pX1Z0) * pX1Z0 * (1 - pZt0)
  Jz <- (1 - pX1Z1) * pX1Z1 * pZt0 - (1 - pX0Z1) * pX0Z1 * pZt0
  J.omtm1.noX <- matrix(c(J0, Jt, Jxt, Jz), nrow = 1)
  RDest.omtm1.noX <- (pX1Z1 * pZt0 - pX0Z1 * pZt0) +
    (pX1Z0 * (1 - pZt0) - pX0Z0 * (1 - pZt0))
  RDse.omtm1.noX <- c(sqrt(J.omtm1.noX %*% vcov.omtm1.noX %*% t(J.omtm1.noX)))
  cat(paste("OMTM1-C convergence code =", code.omtm1.noX, "\n"))
  
  ##### OMTM2 #####
  cat("Fitting OMTM2\n")
  XMat.ppo2 <- XMat.dat[, c(2, 2, 3)]
  # Flip the outcome to fit OMTM2
  Yomtm2.dat <- abs(Y.dat - 5)
  params.omtm2 <- c(-rev(alpha), -beta[1], -(beta[2] + beta.ppo[2]), -beta[3], -beta[4],
                    beta.ppo[2] - beta.ppo[1], beta.ppo[2], 0,
                    rev(c(t(gamma.mat[[1]]))))
  stepmax <- 1
  code <- 5
  while(code == 5){
    cat(paste("Try stepmax =", stepmax, "\n"))
    reset_step <- FALSE
    tryCatch(mdl.omtm2 <- omtm1(
      params = params.omtm2, yval = Yomtm2.dat, XMat = XMat.dat, id = id.dat,
      ppo.k = c(2, 3, 3), x.ppo = XMat.ppo2, u = NULL, ref.muc = 2,
      TransIndMtx = NA, ProfileCol = NA, print.level = 1,
      stepmax = stepmax, UseGrad = TRUE, iterlim = 250, check.analyticals = FALSE),
      error = function(e){reset_step <<- TRUE})
    cat(paste("Reset step", reset_step, "\n"))
    if(reset_step == TRUE){
      stepmax <- stepmax * 0.25
    } else {
      stepmax <- stepmax * 1.5
      params.omtm2 <- c(mdl.omtm2$alpha, mdl.omtm2$beta, c(t(mdl.omtm2$gamma[[1]])))
      code <- code.omtm2 <- mdl.omtm2$control["convergence_code"]
    }}
  coefs.omtm2 <- c(mdl.omtm2$alpha, mdl.omtm2$beta, c(t(mdl.omtm2$gamma[[1]])))
  names(coefs.omtm2) <- c("est.a1", "est.a2", "est.a3",
                          "est.bx1", "est.bt1", "est.bxt1", "est.bz",
                          "est.bt2", "est.bt3", "est.bxt2",
                          "est.g11", "est.g12", "est.g13",
                          "est.g21", "est.g22", "est.g23",
                          "est.g31", "est.g32", "est.g33")
  vcov.omtm2 <- mdl.omtm2$vcov
  ses.omtm2 <- c(sqrt(diag(vcov.omtm2)),
                 vcov.omtm2[4, 6], vcov.omtm2[4, 10], vcov.omtm2[6, 10])
  names(ses.omtm2) <- c("se.a1", "se.a2", "se.a3",
                        "se.bx1", "se.bt1", "se.bxt1", "se.bz",
                        "se.bt2", "se.bt3", "se.bxt2",
                        "se.g12", "se.g13",
                        "se.g22", "se.g23",
                        "se.g32", "se.g33",
                        "cov.bx1bxt1", "cov.bx1bxt2", "cov.bxt1bxt2")
  vcov.omtm2 <- vcov.omtm2[c(3, 4:7, 9:10), c(3, 4:7, 9:10)]
  # Calculate probabilities
  pX1Z1 <- expit(-sum(coefs.omtm2[paste0("est.", c("a3", "bx1", "bt1", "bxt1", "bt3", "bxt2", "bz"))]))
  pX1Z0 <- expit(-sum(coefs.omtm2[paste0("est.", c("a3", "bx1", "bt1", "bxt1", "bt3", "bxt2"))]))
  pX0Z1 <- expit(-sum(coefs.omtm2[paste0("est.", c("a3", "bt1", "bt3", "bz"))]))
  pX0Z0 <- expit(-sum(coefs.omtm2[paste0("est.", c("a3", "bt1", "bt3"))]))
  # Calculate the Jacobian matrix
  J0 <- Jt1 <- Jt3 <- ((1 - pX1Z1) * pX1Z1 * pZt0 - (1 - pX0Z1) * pX0Z1 * pZt0) +
    ((1 - pX1Z0) * pX1Z0 * (1 - pZt0) - (1 - pX0Z0) * pX0Z0 * (1 - pZt0))
  Jx <- Jxt1 <- Jxt2 <- (1 - pX1Z1) * pX1Z1 * pZt0 + (1 - pX1Z0) * pX1Z0 * (1 - pZt0)
  Jz <- (1 - pX1Z1) * pX1Z1 * pZt0 - (1 - pX0Z1) * pX0Z1 * pZt0
  J.omtm2 <- matrix(c(J0, Jx, Jt1, Jxt1, Jz, Jt3, Jxt2), nrow = 1)
  RDest.omtm2 <- (pX1Z1 * pZt0 - pX0Z1 * pZt0) +
    (pX1Z0 * (1 - pZt0) - pX0Z0 * (1 - pZt0))
  RDse.omtm2 <- c(sqrt(J.omtm2 %*% vcov.omtm2 %*% t(J.omtm2)))
  cat(paste("OMTM2 convergence code =", code.omtm2, "\n"))
  
  ##### OMTM2-C #####
  cat("Fitting OMTM2-C \n")
  XMat.ppo2.noX <- XMat.dat.noX[, c(1, 1, 2)]
  params.omtm2.noX <- c(-rev(alpha), -(beta[2] + beta.ppo[2]), -beta[3], -beta[4],
                        beta.ppo[2] - beta.ppo[1], beta.ppo[2], 0,
                        rev(c(t(gamma.mat[[1]]))))
  stepmax <- 1
  code <- 5
  while(code == 5){
    cat(paste("Try stepmax =", stepmax, "\n"))
    reset_step <- FALSE
    tryCatch(
      mdl.omtm2.noX <- omtm1(
        params = params.omtm2.noX, yval = Yomtm2.dat, XMat = XMat.dat.noX, id = id.dat,
        ppo.k = c(2, 3, 3), x.ppo = XMat.ppo2.noX, u = NULL, ref.muc = 2,
        TransIndMtx = NA, ProfileCol = NA, print.level = 1,
        stepmax = stepmax, UseGrad = TRUE, iterlim = 250, check.analyticals = FALSE),
      error = function(e){reset_step <<- TRUE})
    cat(paste("Reset step", reset_step, "\n"))
    if(reset_step == TRUE){
      stepmax <- stepmax * 0.25
    } else {
      stepmax <- stepmax * 1.5
      params.omtm2.noX <- c(mdl.omtm2.noX$alpha, mdl.omtm2.noX$beta, c(t(mdl.omtm2.noX$gamma[[1]])))
      code <- code.omtm2.noX <- mdl.omtm2.noX$control["convergence_code"]
    }}
  coefs.omtm2.noX <- c(mdl.omtm2.noX$alpha, mdl.omtm2.noX$beta, c(t(mdl.omtm2.noX$gamma[[1]])))
  names(coefs.omtm2.noX) <- c("est.a1", "est.a2", "est.a3",
                              "est.bt1", "est.bxt1", "est.bz",
                              "est.bt2", "est.bt3", "est.bxt2",
                              "est.g11", "est.g12", "est.g13",
                              "est.g21", "est.g22", "est.g23",
                              "est.g31", "est.g32", "est.g33")
  vcov.omtm2.noX <- mdl.omtm2.noX$vcov
  ses.omtm2.noX <- c(sqrt(diag(vcov.omtm2.noX)), vcov.omtm2.noX[5, 9])
  names(ses.omtm2.noX) <- c("se.a1", "se.a2", "se.a3",
                            "se.bt1", "se.bxt1", "se.bz",
                            "se.bt2", "se.bt3", "se.bxt2",
                            "se.g12", "se.g13",
                            "se.g22", "se.g23",
                            "se.g32", "se.g33", "cov.bxt1bxt2")
  vcov.omtm2.noX <- vcov.omtm2.noX[c(3, 4:6, 8:9), c(3, 4:6, 8:9)]
  # Calculate probabilities
  pX1Z1 <- expit(-sum(coefs.omtm2.noX[paste0("est.", c("a3", "bt1", "bxt1", "bt3", "bxt2", "bz"))]))
  pX1Z0 <- expit(-sum(coefs.omtm2.noX[paste0("est.", c("a3", "bt1", "bxt1", "bt3", "bxt2"))]))
  pX0Z1 <- expit(-sum(coefs.omtm2.noX[paste0("est.", c("a3", "bt1", "bt3", "bz"))]))
  pX0Z0 <- expit(-sum(coefs.omtm2.noX[paste0("est.", c("a3", "bt1", "bt3"))]))
  # Calculate the Jacobian matrix
  J0 <- Jt1 <- Jt3 <- ((1 - pX1Z1) * pX1Z1 * pZt0 - (1 - pX0Z1) * pX0Z1 * pZt0) +
    ((1 - pX1Z0) * pX1Z0 * (1 - pZt0) - (1 - pX0Z0) * pX0Z0 * (1 - pZt0))
  Jxt1 <- Jxt2 <- (1 - pX1Z1) * pX1Z1 * pZt0 + (1 - pX1Z0) * pX1Z0 * (1 - pZt0)
  Jz <- (1 - pX1Z1) * pX1Z1 * pZt0 - (1 - pX0Z1) * pX0Z1 * pZt0
  J.omtm2.noX <- matrix(c(J0, Jt1, Jxt1, Jz, Jt3, Jxt2), nrow = 1)
  RDest.omtm2.noX <- (pX1Z1 * pZt0 - pX0Z1 * pZt0) +
    (pX1Z0 * (1 - pZt0) - pX0Z0 * (1 - pZt0))
  RDse.omtm2.noX <- c(sqrt(J.omtm2.noX %*% vcov.omtm2.noX %*% t(J.omtm2.noX)))
  cat(paste("OMTM2-C convergence code =", code.omtm2, "\n"))
  
  ##### OMTM3 #####
  cat("Fitting OMTM3\n")
  XMat.ppo3 <- XMat.dat[, c(2, 2, 3, 3)]
  params.omtm3 <- c(alpha, beta, beta.ppo[1:2], 0, 0, c(t(gamma.mat[[1]])))
  stepmax <- 1
  code <- 5
  while(code == 5){
    cat(paste("Try stepmax =", stepmax, "\n"))
    reset_step <- FALSE
    tryCatch(
      mdl.omtm3 <- omtm1(
        params = params.omtm3, yval = Y.dat, XMat = XMat.dat, id = id.dat,
        ppo.k = c(2, 3, 2, 3), x.ppo = XMat.ppo3, u = NULL, ref.muc = NA,
        TransIndMtx = NA, ProfileCol = NA, print.level = 1,
        stepmax = stepmax, UseGrad = TRUE, iterlim = 250, check.analyticals = FALSE),
      error = function(e){reset_step <<- TRUE})
    cat(paste("Reset step", reset_step, "\n"))
    if(reset_step == TRUE){
      stepmax <- stepmax * 0.25
    } else {
      stepmax <- stepmax * 1.5
      params.omtm3 <- c(mdl.omtm3$alpha, mdl.omtm3$beta, c(t(mdl.omtm3$gamma[[1]])))
      code <- code.omtm3 <- mdl.omtm3$control["convergence_code"]
    }}
  coefs.omtm3 <- c(mdl.omtm3$alpha, mdl.omtm3$beta, c(t(mdl.omtm3$gamma[[1]])))
  names(coefs.omtm3) <- c("est.a1", "est.a2", "est.a3",
                          "est.bx1", "est.bt1", "est.bxt1", "est.bz",
                          "est.bt2", "est.bt3", "est.bxt2", "est.bxt3",
                          "est.g11", "est.g12", "est.g13",
                          "est.g21", "est.g22", "est.g23",
                          "est.g31", "est.g32", "est.g33")
  vcov.omtm3 <- mdl.omtm3$vcov
  ses.omtm3 <- c(sqrt(diag(vcov.omtm3)), vcov.omtm3[4, 6])
  names(ses.omtm3) <- c("se.a1", "se.a2", "se.a3",
                        "se.bx1", "se.bt1", "se.bxt1", "se.bz",
                        "se.bt2", "se.bt3", "se.bxt2", "se.bxt3",
                        "se.g12", "se.g13",
                        "se.g22", "se.g23",
                        "se.g32", "se.g33", "cov.bx1bxt1")
  vcov.omtm3 <- vcov.omtm3[c(1, 4:7), c(1, 4:7)]
  # Calculate probabilities
  pX1Z1 <- expit(sum(coefs.omtm3[paste0("est.", c("a1", "bx1", "bt1", "bxt1", "bz"))]))
  pX1Z0 <- expit(sum(coefs.omtm3[paste0("est.", c("a1", "bx1", "bt1", "bxt1"))]))
  pX0Z1 <- expit(sum(coefs.omtm3[paste0("est.", c("a1", "bt1", "bz"))]))
  pX0Z0 <- expit(sum(coefs.omtm3[paste0("est.", c("a1", "bt1"))]))
  # Calculate the Jacobian matrix
  J0 <- Jt <- ((1 - pX1Z1) * pX1Z1 * pZt0 - (1 - pX0Z1) * pX0Z1 * pZt0) +
    ((1 - pX1Z0) * pX1Z0 * (1 - pZt0) - (1 - pX0Z0) * pX0Z0 * (1 - pZt0))
  Jx <- Jxt <- (1 - pX1Z1) * pX1Z1 * pZt0 + (1 - pX1Z0) * pX1Z0 * (1 - pZt0)
  Jz <- (1 - pX1Z1) * pX1Z1 * pZt0 - (1 - pX0Z1) * pX0Z1 * pZt0
  J.omtm3 <- matrix(c(J0, Jx, Jt, Jxt, Jz), nrow = 1)
  RDest.omtm3 <- (pX1Z1 * pZt0 - pX0Z1 * pZt0) +
    (pX1Z0 * (1 - pZt0) - pX0Z0 * (1 - pZt0))
  RDse.omtm3 <- c(sqrt(J.omtm3 %*% vcov.omtm3 %*% t(J.omtm3)))
  cat(paste("OMTM3 convergence code =", code.omtm3, "\n"))
  
  ##### OMTM3-C #####
  cat("Fitting OMTM3-C\n")
  XMat.ppo3.noX <- XMat.dat.noX[, c(1, 1, 2, 2)]
  params.omtm3.noX <- c(alpha, beta[-1], beta.ppo[1:2], 0, 0, c(t(gamma.mat[[1]])))
  stepmax <- 1
  code <- 5
  while(code == 5){
    cat(paste("Try stepmax =", stepmax, "\n"))
    reset_step <- FALSE
    tryCatch(
      mdl.omtm3.noX <- omtm1(
        params = params.omtm3.noX, yval = Y.dat, XMat = XMat.dat.noX, id = id.dat,
        ppo.k = c(2, 3, 2, 3), x.ppo = XMat.ppo3.noX, u = NULL, ref.muc = NA,
        TransIndMtx = NA, ProfileCol = NA, print.level = 1,
        stepmax = stepmax, UseGrad = TRUE, iterlim = 250, check.analyticals = FALSE),
      error = function(e){reset_step <<- TRUE})
    cat(paste("Reset step", reset_step, "\n"))
    if(reset_step == TRUE){
      stepmax <- stepmax * 0.25
    } else {
      stepmax <- stepmax * 1.5
      params.omtm3.noX <- c(mdl.omtm3.noX$alpha, mdl.omtm3.noX$beta, c(t(mdl.omtm3.noX$gamma[[1]])))
      code <- code.omtm3.noX <- mdl.omtm3.noX$control["convergence_code"]
    }}
  coefs.omtm3.noX <- c(mdl.omtm3.noX$alpha, mdl.omtm3.noX$beta, c(t(mdl.omtm3.noX$gamma[[1]])))
  names(coefs.omtm3.noX) <- c("est.a1", "est.a2", "est.a3",
                              "est.bt1", "est.bxt1", "est.bz",
                              "est.bt2", "est.bt3", "est.bxt2", "est.bxt3",
                              "est.g11", "est.g12", "est.g13",
                              "est.g21", "est.g22", "est.g23",
                              "est.g31", "est.g32", "est.g33")
  vcov.omtm3.noX <- mdl.omtm3.noX$vcov
  ses.omtm3.noX <- sqrt(diag(vcov.omtm3.noX))
  names(ses.omtm3.noX) <- c("se.a1", "se.a2", "se.a3",
                            "se.bt1", "se.bxt1", "se.bz",
                            "se.bt2", "se.bt3", "se.bxt2", "se.bxt3",
                            "se.g12", "se.g13",
                            "se.g22", "se.g23",
                            "se.g32", "se.g33")
  vcov.omtm3.noX <- vcov.omtm3.noX[c(1, 4:6), c(1, 4:6)]
  # Calculate probabilities
  pX1Z1 <- expit(sum(coefs.omtm3.noX[paste0("est.", c("a1", "bt1", "bxt1", "bz"))]))
  pX1Z0 <- expit(sum(coefs.omtm3.noX[paste0("est.", c("a1", "bt1", "bxt1"))]))
  pX0Z1 <- expit(sum(coefs.omtm3.noX[paste0("est.", c("a1", "bt1", "bz"))]))
  pX0Z0 <- expit(sum(coefs.omtm3.noX[paste0("est.", c("a1", "bt1"))]))
  # Calculate the Jacobian matrix
  J0 <- Jt <- ((1 - pX1Z1) * pX1Z1 * pZt0 - (1 - pX0Z1) * pX0Z1 * pZt0) +
    ((1 - pX1Z0) * pX1Z0 * (1 - pZt0) - (1 - pX0Z0) * pX0Z0 * (1 - pZt0))
  Jxt <- (1 - pX1Z1) * pX1Z1 * pZt0 + (1 - pX1Z0) * pX1Z0 * (1 - pZt0)
  Jz <- (1 - pX1Z1) * pX1Z1 * pZt0 - (1 - pX0Z1) * pX0Z1 * pZt0
  J.omtm3.noX <- matrix(c(J0, Jt, Jxt, Jz), nrow = 1)
  RDest.omtm3.noX <- (pX1Z1 * pZt0 - pX0Z1 * pZt0) +
    (pX1Z0 * (1 - pZt0) - pX0Z0 * (1 - pZt0))
  RDse.omtm3.noX <- c(sqrt(J.omtm3.noX %*% vcov.omtm3.noX %*% t(J.omtm3.noX)))
  cat(paste("OMTM3-C convergence code =", code.omtm1.noX, "\n"))
  
  ##### Survival ####
  cat("Fitting Cox PH\n")
  # Fit Cox PH model and calculate exp(eta)
  mdl.cox <- coxph(Surv(time, Ybin) ~ X + Z, data = dat.cox, ties = "breslow")
  new.dat <- expand.grid(X = c(0, 1), Z = c(0, 1))
  surv.breslow <- survfit(mdl.cox, newdata = new.dat)
  beta_hat <- coef(mdl.cox)
  dat.cox$exp_eta <- exp(c(t(beta_hat) %*% t(dat.cox[, c("X", "Z")])))
  lambda0_hat <- sapply(0:27, function(t) {
    dat.tmp1 <- dat.cox %>% filter(time == t & Ybin == 1)
    dat.tmp0 <- dat.cox %>% filter(time >= t)
    nrow(dat.tmp1) / sum(dat.tmp0$exp_eta)
  })
  # Breslow estimator for cumulative baseline hazard
  Lambda0_hat <- cumsum(lambda0_hat)
  # Variance
  # for beta_hat
  beta_hat.var <- vcov(mdl.cox)
  # for Lambda0_hat
  d <- c(); tmp0 <- c(); tmp1x <- c(); tmp1z <- c()
  for (t in 0:27) {
    dat.tmp1 <- dat.cox %>% filter(time == t & Ybin == 1)
    dat.tmp0 <- dat.cox %>% filter(time >= t)
    tmp0 <- c(tmp0, sum(dat.tmp0$exp_eta))
    tmp1x <- c(tmp1x, sum(dat.tmp0$X * dat.tmp0$exp_eta))
    tmp1z <- c(tmp1z, sum(dat.tmp0$Z * dat.tmp0$exp_eta))
    d <- c(d, nrow(dat.tmp1))
  }
  Lambda0_hat.var <- c(); cov_hat <- matrix(nrow = mi, ncol = 2)
  for (t in 0:27) {
    mat.tmp <- matrix(c(sum(tmp1x[1:(t+1)] / tmp0[1:(t+1)]^2 * d[1:(t+1)]),
                        sum(tmp1z[1:(t+1)] / tmp0[1:(t+1)]^2 * d[1:(t+1)])), nrow = 1)
    Lambda0_hat.var <- c(Lambda0_hat.var,
                         c(mat.tmp %*% beta_hat.var %*% t(mat.tmp)) + sum(d[1:(t+1)] / tmp0[1:(t+1)]^2))
    # for Cov(Lambda0_hat, beta_hat)
    cov_hat[t+1, ] <- - mat.tmp %*% beta_hat.var
  }
  # Calculate cumulative hazard function at t
  t <- 27
  LambdaX0Z0 <- Lambda0_hat[t+1]
  LambdaX1Z1 <- LambdaX0Z0 * exp(sum(beta_hat))
  LambdaX1Z0 <- LambdaX0Z0 * exp(beta_hat[1])
  LambdaX0Z1 <- LambdaX0Z0 * exp(beta_hat[2])
  # Calculate survival function at t
  SX0Z0 <- exp(-LambdaX0Z0)
  SX1Z1 <- exp(-LambdaX1Z1)
  SX1Z0 <- exp(-LambdaX1Z0)
  SX0Z1 <- exp(-LambdaX0Z1)
  # RD estimation
  RDest.surv <- (SX0Z1 * pZt0 - SX1Z1 * pZt0) + (SX0Z0 * (1 - pZt0) - SX1Z0 * (1 - pZt0))
  # Covariance matrix
  Sigmat_hat <- rbind(matrix(c(Lambda0_hat.var[t+1], cov_hat[t+1,] ), nrow = 1),
                      cbind(t(cov_hat[t+1, , drop = F]), beta_hat.var))
  # Jacobian matrix
  Jt <- matrix(c(
    (-exp(beta_hat[2]) * SX0Z1 * pZt0 + exp(sum(beta_hat)) * SX1Z1 * pZt0) +
      (-SX0Z0 * (1 - pZt0) + exp(beta_hat[1]) * SX1Z0 * (1 - pZt0)),
    LambdaX0Z0 * (exp(sum(beta_hat)) * SX1Z1 * pZt0 + exp(beta_hat[1]) * SX1Z0 * (1 - pZt0)),
    LambdaX0Z0 * (exp(sum(beta_hat)) * SX1Z1 * pZt0 - exp(beta_hat[2]) * SX0Z1 * pZt0)), nrow = 1)
  RDse.surv <- c(sqrt(Jt %*% Sigmat_hat %*% t(Jt)))
  
  ##### Save #####
  cat("Saving Results\n")
  # log(OR)
  coefs.rst <- bind_rows(
    coefs.glmb, coefs.bmtm, coefs.bmtm.noX, coefs.glmo1, coefs.glmo2, coefs.glmo3,
    coefs.omtm1, coefs.omtm1.noX, coefs.omtm2, coefs.omtm2.noX, coefs.omtm3, coefs.omtm3.noX)
  ses.rst <- bind_rows(
    ses.glmb, ses.bmtm, ses.bmtm.noX, ses.glmo1, ses.glmo2, ses.glmo3,
    ses.omtm1, ses.omtm1.noX, ses.omtm2, ses.omtm2.noX, ses.omtm3, ses.omtm3.noX)
  sim1rst <- data.frame(
    data_gen = "PPO T and Z",
    fu       = "Complete",
    mdl_fit  = c("GLMB", "BMTM", "BMTM-C", "GLMO1", "GLMO2", "GLMO3",
                 "OMTM1", "OMTM1-C", "OMTM2", "OMTM2-C", "OMTM3", "OMTM3-C")
  )
  sim1rst <- cbind(sim1rst, coefs.rst, ses.rst)
  if (file.exists("output-missZ")){
    save(sim1rst, file = paste0("output-missZ/sim1-run-", run, ".RData"))}
  # RD
  sim2rst <- data.frame(
    data_gen = "PPO T and Z",
    fu       = "Complete",
    mdl_fit  = c("GLMB", "BMTM", "BMTM-C", "GLMO1", "GLMO2", "GLMO3",
                 "OMTM1", "OMTM1-C", "OMTM2", "OMTM2-C", "OMTM3", "OMTM3-C", "Survival"),
    RDest    = c(RDest.glmb, RDest.bmtm, RDest.bmtm.noX, RDest.glmo1, RDest.glmo2, RDest.glmo3,
                 RDest.omtm1, RDest.omtm1.noX, RDest.omtm2, RDest.omtm2.noX, RDest.omtm3, RDest.omtm3.noX, RDest.surv),
    RDse     = c(RDse.glmb, RDse.bmtm, RDse.bmtm.noX, RDse.glmo1, RDse.glmo2, RDse.glmo3,
                 RDse.omtm1, RDse.omtm1.noX, RDse.omtm2, RDse.omtm2.noX, RDse.omtm3, RDse.omtm3.noX, RDse.surv),
    code     = c(NA, code.bmtm, code.bmtm.noX, rep(NA, 3),
                 code.omtm1, code.omtm1.noX, code.omtm2, code.omtm2.noX, code.omtm3, code.omtm3.noX, NA),
    pZt0     = pZt0
  )
  if (file.exists("output-missZ")){
    save(sim2rst, file = paste0("output-missZ/sim2-run-", run, ".RData"))}
}
