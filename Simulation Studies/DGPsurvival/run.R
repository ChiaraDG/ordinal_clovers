##### This file executes a single simulation and saves data and results
# Load packages
library(OMTM1)  
library(dplyr)
library(survival)
# Check the version of loaded packages
print(packageVersion("OMTM1"));
print(packageVersion("dplyr"));
print(packageVersion("survival"))
# Load simulation design setup
source("setup.R")

simulation <- function(run, data_gen){
  
  ##### Data - Complete ####
  cat("Generating random data from known parameters\n")
  # Generate binary treatment indicator X (`tx`)
  tx <- rbinom(N, size = 1, prob = pX)
  # Generate binary baseline status Z (Fix p(Z) = 0.5)
  Z <- rep(0, N)
  Z[sample(1:N, N * pZ)] <- 1
  if(data_gen == "Survival"){
    # Construct the design matrix `XMat` in the form of (tx, Z)
    XMat <- cbind(tx, Z)
    # Generate survival time 
    lp <- c(beta %*% t(XMat))
    time <- ceiling((-(log(runif(N)) / (lambda*exp(lp))))^(1/nu))
    Ybin <- ifelse(time <= max(t), 1, 0)
    time <- ifelse(time > max(t), max(t), time)
    id <- 1:N
    # Incorporate the generated data into a dataframe
    dat.cox <- data.frame(id = id, X = tx, time = time, Z = Z, Ybin = Ybin)
    # Create a dataframe for BMTM
    dat <- data.frame()
    # Iterate over each row in the original dataframe
    for (i in 1:nrow(dat.cox)) {
      # Extract the current row
      tmpi <- dat.cox[i, ]
      # Create a temp dataframe with mi rows for the current id
      tmp_df <- data.frame(
        id = rep(tmpi$id, mi),
        X = rep(tmpi$X, mi),
        time = 0:(mi - 1),
        Z = rep(tmpi$Z, mi),
        Ybin = rep(2, mi)  # Initialize Ybin as 2
      ) %>%
        mutate(logt = log(time + 1, base = mi))
      # If Ybin in the original dataframe is 1, update Ybin in tmp_df accordingly
      if (tmpi$Ybin == 1) {
        tmp_df$Ybin[ceiling(tmpi$time + 1):mi] <- 1
      }
      # Append the temporary dataframe to the result dataframe
      dat <- rbind(dat, tmp_df)
    }
    dat0 <- dat
  }
  # Save data
  if (file.exists("data-complete")){
    save(dat.cox, file = paste0("data-complete/data-", run, ".RData"))}
  
  ## Calculate p(Z) at the beginning of the study
  dat.t0 <- dat %>% filter(logt == 0)
  pZt0 <- mean(dat.t0$Z == 1)

  ##### GLMB ####
  cat("Fitting GLMB\n")
  # Subset data at the end of follow-up time only
  dat.cross <- subset(dat, logt == 1)
  # Create binary outcome variable: Ybin=1 (event Y=1), Ybin=0 (no event Y=2:4)
  dat.cross$Ybin <- ifelse(dat.cross$Y == 1 & !is.na(dat.cross$Y), 1, 0)
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
  # Fit BMTM model and extract coefficients, covariance matrix
  XMat.dat <- with(dat, cbind(X, logt, X*logt, Z))
  Ybin.dat <- dat$Ybin
  id.dat <- dat$id
  params.bmtm <- c(-4.75, 0, 2, beta, 25)
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
  params.bmtm.noX <- c(-4.75, 2, beta, 25)
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

  ##### Survival ####
  cat("Fitting Cox PH\n")
  # Fit Cox PH model and calculate exp(eta)
  mdl.cox <- coxph(Surv(time, Ybin) ~ X + Z, data = dat.cox, ties = "breslow")
  new.dat <- expand.grid(X = c(0, 1), Z = c(0, 1))
  surv.breslow <- survfit(mdl.cox, newdata = new.dat)
  beta_hat <- coef(mdl.cox)
  dat.cox$exp_eta <- exp(c(t(beta_hat) %*% t(dat.cox[, c("X", "Z")])))
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
    cov_hat[t+1,] <- - mat.tmp %*% beta_hat.var
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
  # RD
  sim2rst <- data.frame(
    data_gen = data_gen,
    fu       = "Complete",
    mdl_fit  = c("GLMB", "BMTM", "BMTM-C", "Survival"),
    RDest    = c(RDest.glmb, RDest.bmtm, RDest.bmtm.noX, RDest.surv),
    RDse     = c(RDse.glmb, RDse.bmtm, RDse.bmtm.noX, RDse.surv),
    code     = c(NA, code.bmtm, code.bmtm.noX, NA),
    pZt0     = pZt0
  )
  if (file.exists("output-complete")){
    save(sim2rst, file = paste0("output-complete/sim2-run-", run, ".RData"))}

  #### Data - DAR (mild) ####
  dat <- dat0
  prob.miss <- 0.01
  # Use a list to store each subject's data
  dat.missing <- vector("list", N)  
  for (i in 1:N) {  # For each subject i
    id.dat <- subset(dat, id == i)
    Y.missing <- id.dat$Ybin  # Start with the full outcome vector
    ## Generate dropout indicators
    missing.inds <- rbinom(n = mi - 1, size = 1, prob = prob.miss)
    missing.inds <- c(0, missing.inds)  # Ensure the first outcome is always observed
    for (j in 2:mi) {
      if (is.na(Y.missing[j - 1])) {  # If the previous outcome was missing
        Y.missing[j] <- NA
      } else if (missing.inds[j] == 1) {  # If dropout occurs
        Y.missing[j] <- NA
      }
      if (is.na(Y.missing[j]) && !is.na(Y.missing[j - 1]) && Y.missing[j - 1] == 1) {
        Y.missing[j] <- 1  # Set to absorbing state if applicable (even if dropout occurs)
      }
    }
    dat.missing[[i]] <- data.frame(id.dat, Y.missing = Y.missing)
  }
  # Combine all subjects' data into a single data frame
  dat.missing <- do.call(rbind, dat.missing)
  # Filter out rows with missing outcomes and remove the Y.missing column
  dat <- dat.missing %>% filter(!is.na(Y.missing)) %>% select(-Y.missing)
  dat.tmp <- dat
  # Create binary outcome variable: Ybin=1 (event Y=1), Ybin=0 (no event Y=2:4)
  dat.tmp$Ybin <- ifelse(dat.tmp$Ybin == 1, 1, 0)
  id.event <- dat.tmp[with(dat.tmp, logt == 1 & Ybin == 1),]$id
  dat.noevent <- dat.tmp %>% filter(!(id %in% id.event)) %>%
    group_by(id) %>% filter(time == max(time)) %>% ungroup()
  dat.event <- dat.tmp %>% filter(id %in% id.event & Ybin == 1) %>%
    group_by(id) %>% filter(time == min(time)) %>% ungroup()
  dat.cox <- rbind(dat.noevent, dat.event)
  if (file.exists("data-dar-mild")){
    save(dat.cox, file = paste0("data-dar-mild/data-", run, ".RData"))}
  
  # Calculate p(Z) at the beginning of the study 
  dat.t0 <- dat %>% filter(logt == 0)
  pZt0 <- mean(dat.t0$Z == 1)
  
  ##### GLMB ####
  cat("Fitting GLMB\n")
  # Subset data at the end of follow-up time only
  dat.cross <- subset(dat, logt == 1)
  # Create binary outcome variable: Ybin=1 (event Y=1), Ybin=0 (no event Y=2:4)
  dat.cross$Ybin <- ifelse(dat.cross$Y == 1 & !is.na(dat.cross$Y), 1, 0)
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
  # Fit BMTM model and extract coefficients, covariance matrix
  XMat.dat <- with(dat, cbind(X, logt, X*logt, Z))
  Ybin.dat <- dat$Ybin
  id.dat <- dat$id
  params.bmtm <- c(-4.75, 0, 2, beta, 25)
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
  params.bmtm.noX <- c(-4.75, 2, beta, 25)
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
  sim2rst <- data.frame(
    data_gen = data_gen,
    fu       = "DAR (mild)",
    mdl_fit  = c("GLMB", "BMTM", "BMTM-C", "Survival"),
    RDest    = c(RDest.glmb, RDest.bmtm, RDest.bmtm.noX, RDest.surv),
    RDse     = c(RDse.glmb, RDse.bmtm, RDse.bmtm.noX, RDse.surv),
    code     = c(NA, code.bmtm, code.bmtm.noX, NA),
    pZt0     = pZt0
  )
  if (file.exists("output-dar-mild")){
    save(sim2rst, file = paste0("output-dar-mild/sim2-run-", run, ".RData"))}
  
  #### Data - DAR (severe) ####
  dat <- dat0
  prob.miss <- 0.04
  # Use a list to store each subject's data
  dat.missing <- vector("list", N)  
  for (i in 1:N) {  # For each subject i
    id.dat <- subset(dat, id == i)
    Y.missing <- id.dat$Ybin  # Start with the full outcome vector
    ## Generate dropout indicators
    missing.inds <- rbinom(n = mi - 1, size = 1, prob = prob.miss)
    missing.inds <- c(0, missing.inds)  # Ensure the first outcome is always observed
    for (j in 2:mi) {
      if (is.na(Y.missing[j - 1])) {  # If the previous outcome was missing
        Y.missing[j] <- NA
      } else if (missing.inds[j] == 1) {  # If dropout occurs
        Y.missing[j] <- NA
      }
      if (is.na(Y.missing[j]) && !is.na(Y.missing[j - 1]) && Y.missing[j - 1] == 1) {
        Y.missing[j] <- 1  # Set to absorbing state if applicable (even if dropout occurs)
      }
    }
    dat.missing[[i]] <- data.frame(id.dat, Y.missing = Y.missing)
  }
  # Combine all subjects' data into a single data frame
  dat.missing <- do.call(rbind, dat.missing)
  # Filter out rows with missing outcomes and remove the Y.missing column
  dat <- dat.missing %>% filter(!is.na(Y.missing)) %>% select(-Y.missing)
  dat.tmp <- dat
  # Create binary outcome variable: Ybin=1 (event Y=1), Ybin=0 (no event Y=2:4)
  dat.tmp$Ybin <- ifelse(dat.tmp$Ybin == 1, 1, 0)
  id.event <- dat.tmp[with(dat.tmp, logt == 1 & Ybin == 1),]$id
  dat.noevent <- dat.tmp %>% filter(!(id %in% id.event)) %>%
    group_by(id) %>% filter(time == max(time)) %>% ungroup()
  dat.event <- dat.tmp %>% filter(id %in% id.event & Ybin == 1) %>%
    group_by(id) %>% filter(time == min(time)) %>% ungroup()
  dat.cox <- rbind(dat.noevent, dat.event)
  if (file.exists("data-dar-severe")){
    save(dat.cox, file = paste0("data-dar-severe/data-", run, ".RData"))}
  
  # Calculate p(Z) at the beginning of the study 
  dat.t0 <- dat %>% filter(logt == 0)
  pZt0 <- mean(dat.t0$Z == 1)
  
  ##### GLMB ####
  cat("Fitting GLMB\n")
  # Subset data at the end of follow-up time only
  dat.cross <- subset(dat, logt == 1)
  # Create binary outcome variable: Ybin=1 (event Y=1), Ybin=0 (no event Y=2:4)
  dat.cross$Ybin <- ifelse(dat.cross$Y == 1 & !is.na(dat.cross$Y), 1, 0)
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
  # Fit BMTM model and extract coefficients, covariance matrix
  XMat.dat <- with(dat, cbind(X, logt, X*logt, Z))
  Ybin.dat <- dat$Ybin
  id.dat <- dat$id
  params.bmtm <- c(-4.75, 0, 2, beta, 25)
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
  params.bmtm.noX <- c(-4.75, 2, beta, 25)
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
  sim2rst <- data.frame(
    data_gen = data_gen,
    fu       = "DAR (severe)",
    mdl_fit  = c("GLMB", "BMTM", "BMTM-C", "Survival"),
    RDest    = c(RDest.glmb, RDest.bmtm, RDest.bmtm.noX, RDest.surv),
    RDse     = c(RDse.glmb, RDse.bmtm, RDse.bmtm.noX, RDse.surv),
    code     = c(NA, code.bmtm, code.bmtm.noX, NA),
    pZt0     = pZt0
  )
  if (file.exists("output-dar-severe")){
    save(sim2rst, file = paste0("output-dar-severe/sim2-run-", run, ".RData"))}
}
