library(tidyverse)
library(OMTM1)
library(VGAM)

##########################################################################################
########## Case1 : One absorbing states log(time) ########################################
##########################################################################################

# load the outcome
dat           <- read.csv("./Data/dailystateICU28daysmortality_oneabssorbing.csv")
dat[dat$ReferenceID == "U01-01301" & dat$Day >=21,"DailyState"] <- "Discharged"
write.csv(dat,"./Data/dailystateICU28daysmortality_oneabssorbing.csv")
# load the covariate
covs          <- read.csv("./Data/baselinecovariates.csv")
# merge the data
dat           <- inner_join(dat, covs, by = "ReferenceID")

# redefine the variables we need
dat           <- subset(dat, Day > 0)
dat$logtime   <- log((dat$Day + 1), base = 29)
dat$trt       <- ifelse(dat$rand_trt == "Restrictive Fluid Group", 1, 0)
dat$sofa      <- (dat$d_sofa_gcs - mean(dat$d_sofa_gcs))/3
dat$age       <- (dat$scr_age - mean(dat$scr_age))/15
dat$female    <- ifelse(dat$dmg_sex == "Female", 1, 0)
dat$ards      <- ifelse(dat$ards == "Yes", 1, 0)
dat$DailyState2 <- ifelse(dat$DailyState == "Death", 1,
                          ifelse(dat$DailyState == "Discharged", 4,
                                 ifelse(dat$DailyState == "Hospital", 3,
                                        2)))

id          <- dat$ReferenceID
Y.0         <- dat$DailyState2


#########################################################################################################
##### Cross-Sectional Outcome ###########################################################################
#########################################################################################################

#########################################################################################################
##### Longitudinal Outcome ##############################################################################
#########################################################################################################
Y.mtlv <- ifelse(dat$DailyState2 == 1, 1, 2)
dat$Y.bin  <- ifelse(dat$DailyState2 == 1, 1, 0)
mod.glm <- glm(Y.bin ~ trt + sofa + age + female,
                data = subset(dat, logtime == 1),
                family="binomial")
summary(mod.glm)

XMat        <- with(dat, cbind(trt, ns(logtime, 3),
                               trt*ns(logtime, 3), sofa, age, female))

XMat.noX    <- with(dat, cbind(ns(logtime, 3),
                               trt*ns(logtime, 3), sofa, age, female))

##### binary case #########################################################################
Y.mtlv <- ifelse(dat$DailyState2 == 1, 1, 2)
dat$Y.bin  <- ifelse(dat$DailyState2 == 1, 1, 0)

# set up starting parameters: linear time
mod.vgam <- glm(Y.bin ~ trt*ns(dat$logtime,3) +
                    sofa + age + female, data = dat,
                family="binomial")

mod.vgam.noX <- glm(Y.bin ~ ns(dat$logtime,3) + trt:ns(dat$logtime,3) +
                    sofa + age + female, data = dat,
                family="binomial")

params        <- c(coef(mod.vgam)[c(1:5, 9:11, 6:8)], 0)
params.noX    <- c(coef(mod.vgam.noX)[1], coef(mod.vgam.noX)[c(2:4, 8:10, 5:7)], 0)
mod <- omtm1(params = params,
             yval = Y.mtlv, XMat = XMat,
             id = id, ppo.k=NULL, x.ppo=NULL, u=NULL,
             ref.muc = 2, TransIndMtx=NA,
             stepmax = 2, UseGrad=TRUE, iterlim=200,
             check.analyticals = FALSE, print.level=2)

mod.noX <- omtm1(params = params.noX,
             yval = Y.mtlv, XMat = XMat.noX,
             id = id, ppo.k=NULL, x.ppo=NULL, u=NULL,
             ref.muc = 2, TransIndMtx=NA,
             stepmax = 2, UseGrad=TRUE, iterlim=200,
             check.analyticals = FALSE, print.level=2)

saveRDS(mod, file="bmtm.RData")
saveRDS(mod.noX, file="bmtm_noX.RData")

##### ordinal case #########################################################################
##### PPO for time and time by treatment #####
XMat.ppo    <- with(dat, cbind(ns(logtime, 3),
                               ns(logtime, 3),
                               trt*ns(logtime, 3), trt*ns(logtime,3)))

mod.vgam <- vglm(DailyState2 ~ trt*ns(logtime, 3) + sofa + age + female, data = dat,
                 family=cumulative(parallel=TRUE,
                                   reverse = FALSE))

mod.vgam.noX <- vglm(DailyState2 ~ ns(logtime, 3) + trt:ns(logtime, 3) + sofa + age + female, data = dat,
                 family=cumulative(parallel=TRUE,
                                   reverse = FALSE))
gamma.mat <- list()
gamma.mat[[1]] <- rbind(c(0, 0, 0),
                        c(0, 0, 0),
                        c(0, 0, 0))
params    <- c(coef(mod.vgam)[c(1:7, 11:13, 8:10)], rep(0, 12),
               c(t(gamma.mat[[1]])))#, t(gamma.mat[[2]])))
params.noX    <- c(coef(mod.vgam.noX)[c(1:6, 10:12, 7:9)], rep(0, 12),
               c(t(gamma.mat[[1]])))

stepmax <- 1.5
code    <- 5

while(code == 5){

    print(stepmax)

    reset_step <- FALSE

    tryCatch(mod <- omtm1(params = params, yval = Y.0, XMat = XMat,
                          id = id, ppo.k=c(2,2,2,3,3,3,2,2,2,3,3,3), x.ppo=XMat.ppo, u=NULL,
                          ref.muc = 2, TransIndMtx=NA,
                          stepmax = stepmax, UseGrad=TRUE, iterlim=200,
                          check.analyticals = FALSE, print.level=2),
             error = function(e){reset_step <<- TRUE})

    print(reset_step)

    if(reset_step == TRUE){
        #stepmax = 0.001
        stepmax = stepmax*0.25
    } else {
        stepmax <- stepmax*1.5
        params  <- c(mod$alpha, mod$beta, c(t(mod$gamma[[1]])))
        code    <- mod$control["convergence_code"]
        print(code)
    }

}
saveRDS(mod, file="omtm3.RData")

stepmax <- 1.5
code    <- 5
while(code == 5){

    print(stepmax)

    reset_step <- FALSE

    tryCatch(mod <- omtm1(params = params.noX, yval = Y.0, XMat = XMat.noX,
                          id = id, ppo.k=c(2,2,2,3,3,3,2,2,2,3,3,3), x.ppo=XMat.ppo, u=NULL,
                          ref.muc = 2, TransIndMtx=NA,
                          stepmax = stepmax, UseGrad=TRUE, iterlim=200,
                          check.analyticals = FALSE, print.level=2),
             error = function(e){reset_step <<- TRUE})

    print(reset_step)

    if(reset_step == TRUE){
        #stepmax = 0.001
        stepmax = stepmax*0.25
    } else {
        stepmax <- stepmax*1.5
        params.noX  <- c(mod$alpha, mod$beta, c(t(mod$gamma[[1]])))
        code    <- mod$control["convergence_code"]
        print(code)
    }

}

saveRDS(mod, file="omtm3_noX.RData")


##### PPO for time only #####
XMat.ppo    <- with(dat, cbind(ns(logtime, 3),
                               ns(logtime, 3)))
# set up starting parameters: linear time
mod.vgam <- vglm(DailyState2 ~ trt*ns(logtime, 3) + sofa +
                     age + female, data = dat,
                 family=cumulative(parallel=TRUE,
                                   reverse = FALSE))
gamma.mat <- list()
gamma.mat[[1]] <- rbind(c(0, 0, 0),
                        c(0, 0, 0),
                        c(0, 0, 0))
params    <- c(coef(mod.vgam)[c(1:7, 11:13, 8:10)],  rep(0, 6),
               c(t(gamma.mat[[1]])))#, t(gamma.mat[[2]])))
params.noX    <- c(coef(mod.vgam.noX)[c(1:6, 10:12, 7:9)], rep(0, 6),
                   c(t(gamma.mat[[1]])))



stepmax <- 1.5
code    <- 5

while(code == 5){

    print(stepmax)

    reset_step <- FALSE

    tryCatch(mod <- omtm1(params = params, yval = Y.0, XMat = XMat,
                          id = id, ppo.k=c(2,2,2,3,3,3), x.ppo=XMat.ppo, u=NULL,
                          ref.muc = 2, TransIndMtx=NA,
                          stepmax = stepmax, UseGrad=TRUE, iterlim=200,
                          check.analyticals = FALSE, print.level=2),
             error = function(e){reset_step <<- TRUE})

    print(reset_step)

    if(reset_step == TRUE){
        #stepmax = 0.001
        stepmax = stepmax*0.25
    } else {
        stepmax <- stepmax*1.5
        params  <- c(mod$alpha, mod$beta, c(t(mod$gamma[[1]])))
        code    <- mod$control["convergence_code"]
        print(code)
    }

}
saveRDS(mod, file="omtm1.RData")

stepmax <- 1.5
code    <- 5

while(code == 5){

    print(stepmax)

    reset_step <- FALSE

    tryCatch(mod <- omtm1(params = params.noX, yval = Y.0, XMat = XMat.noX,
                          id = id, ppo.k=c(2,2,2,3,3,3), x.ppo=XMat.ppo, u=NULL,
                          ref.muc = 2, TransIndMtx=NA,
                          stepmax = stepmax, UseGrad=TRUE, iterlim=200,
                          check.analyticals = FALSE, print.level=2),
             error = function(e){reset_step <<- TRUE})

    print(reset_step)

    if(reset_step == TRUE){
        #stepmax = 0.001
        stepmax = stepmax*0.25
    } else {
        stepmax <- stepmax*1.5
        params.noX  <- c(mod$alpha, mod$beta, c(t(mod$gamma[[1]])))
        code    <- mod$control["convergence_code"]
        print(code)
    }

}
saveRDS(mod, file="omtm1_noX.RData")

#### PPO for time and time by X but only the death category ####
XMat.ppo    <- with(dat, cbind(ns(logtime, 3),
                               ns(logtime, 3),
                               trt*ns(logtime, 3)))

dat$DailyState3 <- ifelse(dat$DailyState2 == 4, 1, NA)
dat$DailyState3 <- ifelse(dat$DailyState2 == 3, 2, dat$DailyState3)
dat$DailyState3 <- ifelse(dat$DailyState2 == 2, 3, dat$DailyState3)
dat$DailyState3 <- ifelse(dat$DailyState2 == 1, 4, dat$DailyState3)
table(dat$DailyState3)
mod.vgam <- vglm(DailyState3 ~ trt*ns(logtime, 3) + sofa +
                     age + female, data = dat,
                 family=cumulative(parallel=TRUE,
                                   reverse = FALSE))
mod.vgam.noX <- vglm(DailyState3 ~ ns(logtime, 3) + trt:ns(logtime, 3) + sofa +
                     age + female, data = dat,
                 family=cumulative(parallel=TRUE,
                                   reverse = FALSE))
gamma.mat <- list()
gamma.mat[[1]] <- rbind(c(0, 0, 0),
                        c(0, 0, 0),
                        c(0, 0, 0))
params    <- c(coef(mod.vgam)[c(1:7, 11:13, 8:10)],  rep(0, 9),
               c(t(gamma.mat[[1]])))#, t(gamma.mat[[2]])))
params.noX    <- c(coef(mod.vgam.noX)[c(1:6, 10:12, 7:9)], rep(0, 9),
                   c(t(gamma.mat[[1]])))

Y.3 <- dat$DailyState3


stepmax <- 1.5
code    <- 5

while(code == 5){

    print(stepmax)

    reset_step <- FALSE

    tryCatch(mod <- omtm1(params = params, yval = Y.3, XMat = XMat,
                          id = id, ppo.k=c(2,2,2,3,3,3,3,3,3), x.ppo=XMat.ppo, u=NULL,
                          ref.muc = 3, TransIndMtx=NA,
                          stepmax = stepmax, UseGrad=TRUE, iterlim=200,
                          check.analyticals = FALSE, print.level=2),
             error = function(e){reset_step <<- TRUE})

    print(reset_step)

    if(reset_step == TRUE){
        #stepmax = 0.001
        stepmax = stepmax*0.25
    } else {
        stepmax <- stepmax*1.5
        params  <- c(mod$alpha, mod$beta, c(t(mod$gamma[[1]])))
        code    <- mod$control["convergence_code"]
        print(code)
    }

}

saveRDS(mod, file="omtm2.RData")

stepmax <- 1.5
code    <- 5
while(code == 5){

    print(stepmax)

    reset_step <- FALSE

    tryCatch(mod <- omtm1(params = params.noX, yval = Y.3, XMat = XMat.noX,
                          id = id, ppo.k=c(2,2,2,3,3,3,3,3,3), x.ppo=XMat.ppo, u=NULL,
                          ref.muc = 3, TransIndMtx=NA,
                          stepmax = stepmax, UseGrad=TRUE, iterlim=200,
                          check.analyticals = FALSE, print.level=2),
             error = function(e){reset_step <<- TRUE})

    print(reset_step)

    if(reset_step == TRUE){
        #stepmax = 0.001
        stepmax = stepmax*0.25
    } else {
        stepmax <- stepmax*1.5
        params.noX  <- c(mod$alpha, mod$beta, c(t(mod$gamma[[1]])))
        code    <- mod$control["convergence_code"]
        print(code)
    }

}
saveRDS(mod, file="omtm2_noX.RData")
