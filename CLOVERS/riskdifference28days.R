##########################################################################################
########## Compute Risk Difference for CLOVERS results ###################################
##########################################################################################

#### load the data we need ###############################################################

library(tidyverse)
library(viridis)
library(VGAM)
library(survival)
library(mvtnorm)
library(msm)

# load the outcome
dat           <- read.csv("/Users/chiaradigravio/Documents/PhD/Collaborations/CLOVERS/Data/dailystateICU28daysmortality_oneabssorbing.csv")
dat[dat$ReferenceID == "U01-01301" & dat$Day >= 21, "DailyState2"] <- "Discharged"
# load the covariate
covs          <- read.csv("/Users/chiaradigravio/Documents/PhD/Collaborations/CLOVERS/Data/baselinecovariates.csv")
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

#### cross-sectional outcomes ###############################################################
dat$Y.glm  <- ifelse(dat$DailyState2 == 1, 1, 0)
dat$Y.mtlv <- ifelse(dat$DailyState2 == 1, 1, 2)
d.cross    <- subset(dat, Day == 28)
# GLM
mod.glm <- glm(Y.glm ~ trt + sofa + age + female, data = d.cross,
               family = "binomial")

dir <- "~/Documents/PhD/Collaborations/CLOVERS/BinaryOrdinal28days/one_absorbing_full/"

# mTLV
mtlv     <- readRDS(paste0(dir, "bmtm.RData"))
mtlv.noX <- readRDS(paste0(dir, "bmtm_noX.RData"))
# VGAM
mod.vgam1 <- vglm(DailyState2 ~ trt + sofa + age + female,
                  data = d.cross,
                  family=cumulative(parallel=TRUE))
mat <- cbind(rbind(1, 0, 0), rbind(0, 1, 1))
clist <- list("(Intercept)" = diag(3), "trt" = mat,
              "sofa" = rbind(1, 1, 1), "age" = rbind(1, 1, 1),
              "female" = rbind(1, 1, 1))
mod.vgam2 <- vglm(DailyState2 ~ trt + sofa + age + female,
                  data = d.cross, constraints = clist,
                  family=cumulative(parallel=FALSE))
clist <- list("(Intercept)" = diag(3), "trt" = diag(3),
              "sofa" = rbind(1, 1, 1), "age" = rbind(1, 1, 1),
              "female" = rbind(1, 1, 1))
mod.vgam3 <- vglm(DailyState2 ~ trt + sofa + age + female,
                  data = d.cross, constraints = clist,
                  family=cumulative(parallel=FALSE))

# OMTM1 and OMTM1-C
omtm1    <- readRDS(paste0(dir, "omtm1.RData"))
omtm1.noX <- readRDS(paste0(dir, "omtm1_noX.RData"))
# OMTM2 and OMTM2-C
omtm2    <- readRDS(paste0(dir, "omtm2.RData"))
omtm2.noX <- readRDS(paste0(dir, "omtm2_noX.RData"))
# OMTM 3 and OMTM3-C
omtm3     <- readRDS(paste0(dir, "omtm3.RData"))
omtm3.noX <- readRDS(paste0(dir, "omtm3_noX.RData"))

# survival model
dat$death_time   <- ifelse(dat$DailyState == "Death", dat$Day, 28)
dat$death_status <-  ifelse(dat$DailyState2 == 1, 1, 0)
dat              <- dat %>% group_by(ReferenceID) %>% mutate(death_time = min(death_time)) %>% ungroup()
dat.surv         <- dat %>% filter(Day == 28)
survs            <- coxph(Surv(death_time, death_status)  ~ trt + sofa + age + female, data = dat.surv,
                          ties = "breslow")

# expit
expit <- function(x){exp(x)/(1 + exp(x))}

#### create a function to compute the risk difference at day 28 for GLM, GLMB and OMTM #########################################################################

rd_computation <- function(coefficient, design_matrix, method){
    # create the design matrix

    d0           <- design_matrix
    d1           <- design_matrix

    if(method %in% c("glm", "vgam1", "survival")){
        # change the treatment effect to 0 and 1
        d0[, "trt"]  <- 0
        d1[, "trt"]  <- 1
    }

    if(method %in% c("vgam2")){
        # change the treatment effect to 0 and 1
        for(i in 1:nrow(design_matrix)){
            d0[, c("trt:1", "trt:2")]  <- 0
            if(d1[i,"(Intercept):1"] == 1){d1[i, "trt:1"] = 1}
            if(d1[i,"(Intercept):2"] == 1){d1[i, "trt:2"] = 1}
            if(d1[i,"(Intercept):3"] == 1){d1[i, "trt:2"] = 1}
        }
    }

    if(method %in% c("vgam3")){
        # change the treatment effect to 0 and 1
        for(i in 1:nrow(design_matrix)){
            d0[, c("trt:1", "trt:2")]  <- 0
            if(d1[i,"(Intercept):1"] == 1){d1[i, "trt:1"] = 1}
            if(d1[i,"(Intercept):2"] == 1){d1[i, "trt:2"] = 1}
            if(d1[i,"(Intercept):3"] == 1){d1[i, "trt:3"] = 1}
        }
    }
    if(method %in% c("bmtm", "omtm1", "omtm3")){
        d0[, c("trt", "trt*t1", "trt*t2", "trt*t3")]  <- 0
        d1[, c("trt")]  <- 1
        d1[, c("trt*t1", "trt*t2", "trt*t3")]  <- d1[, c("t1", "t2", "t3")]
    }

    if(method %in% c("bmtm.noX", "omtm1.noX", "omtm3.noX")){
        d0[, c("trt*t1", "trt*t2", "trt*t3")]  <- 0
        d1[, c("trt*t1", "trt*t2", "trt*t3")]  <- d1[, c("t1", "t2", "t3")]
    }

    if(method == "omtm2"){
        d0[, c("trt", "trt*t1", "trt*t2", "trt*t3", "trt*t1:1", "trt*t2:1", "trt*t3:1")]  <- 0
        d1[, c("trt")]  <- 1
        d1[, c("trt*t1", "trt*t2", "trt*t3", "trt*t1:1", "trt*t2:1", "trt*t3:1")]  <-
            d1[, c("t1", "t2", "t3", "t1", "t2", "t3")]
    }

    if(method == "omtm2.noX"){
        d0[, c("trt*t1", "trt*t2", "trt*t3", "trt*t1:1", "trt*t2:1", "trt*t3:1")]  <- 0
        d1[, c("trt*t1", "trt*t2", "trt*t3", "trt*t1:1", "trt*t2:1", "trt*t3:1")]  <-
            d1[, c("t1", "t2", "t3", "t1", "t2", "t3")]
    }


    # generate new coeffiecient from rmorn
    beta_new <- rmvnorm(1, mean = coefficient$coefs, coefficient$vcovs)

    # compute predicted values
    p1 <- expit(d1%*%t(beta_new))
    p0 <- expit(d0%*%t(beta_new))

    # for coxph
    if(method == "survival"){
        # baseline hazard
        H0 <- basehaz(survs, centered=FALSE)
        p1 <- 1 - exp(-H0[28,1]*exp(d1%*%t(beta_new)))
        p0 <- 1 - exp(-H0[28,1]*exp(d0%*%t(beta_new)))
    }

    if(method %in% c("vgam1", "vgam2", "vgam3")){
        p1 <- matrix(p1, ncol = 3, byrow = TRUE)[,1]
        p0 <- matrix(p0, ncol = 3, byrow = TRUE)[,1]
    }

    rd <- mean(p1 - p0)

    rd

}


#### bootstrap function for survival models: Cox models and multistate models

boot_rd <- function(method, nboot = 10, seed = 123){

    set.seed(seed)
    rd <- c()

    if(method == "coxph"){
        ids   <- unique(dat.surv$ReferenceID)
        for(i in 1:nboot){
            print(i)
            tosample <- sample(1:nrow(dat.surv), size = length(ids), replace = TRUE)
            dat.boot <- dat.surv[tosample,]
            survs    <- coxph(Surv(death_time, death_status)  ~ trt + sofa + age + female, data = dat.boot,
                                      ties = "breslow")
            d1.boot <- dat.boot |> mutate(trt = 1)
            d0.boot <- dat.boot |> mutate(trt = 0)
            bhaz1 <- basehaz(survs, newdata = data.frame(d1.boot), centered = FALSE) |> select(-time)
            bhaz0 <- basehaz(survs, newdata = data.frame(d0.boot), centered = FALSE) |> select(-time)
            p1 <- unlist(bhaz1[28,])*
                exp(as.matrix(d1.boot[,c("trt", "sofa", "age", "female")])%*%survs$coefficients)
            p0 <-  1 - predict(survs, newdata =  d0.boot ,  type = "survival")

            rd[i] <- mean(p1 - p0)
        }
    }

    if(method == "multistate"){
        rd <- c()

        dat.small <- dat |>
            select(ReferenceID, Day, DailyState2, trt, age, sofa, female) |>
            mutate(death_day = ifelse(DailyState2 == 1, 1, 0)) |>
            group_by(ReferenceID) |>
            mutate(death_tot = cumsum(death_day)) |> filter(death_tot <= 1) |>
            mutate(nvisit = n()) |> filter(nvisit > 1) |>
            ungroup() |> data.frame()
        dat.small$obstype <- 2
        Q <- rbind(c(1, 0, 0, 0),
                   c(1, 1, 1, 1),
                   c(1, 1, 1, 1),
                   c(1, 1, 1, 1))
        ids   <- unique(dat.small$ReferenceID)

        for(i in 1:nboot){
            print(i)
            tosample <- sample(ids, size = length(ids), replace = TRUE)
            dat.boot <- dat.small[unlist(lapply(tosample, function(x) which(x == dat.small$ReferenceID))), ]
            # create new ids for the code to run (the msm will not take the duplicated id)
            dat.boot <- dat.boot |> group_by(ReferenceID, Day) |> mutate(nrep = row_number()) |>
                mutate(new_id = paste0(ReferenceID, nrep)) |> ungroup()
            mod.msm <- msm(DailyState2 ~ Day,
                           subject = new_id,
                           data = dat.boot,
                           qmatrix = Q,
                           gen.inits=TRUE,
                           covariates = ~ trt + sofa + female + age,
                           obstype = obstype)
            p1 <- prevalence.msm(mod.msm, covariates=list(trt = 1, female = dat.boot$female,
                                                          sofa = dat.boot$sofa, age = dat.boot$age))$`Expected percentages`[11, "State 1"]/100
            p0 <- prevalence.msm(mod.msm, covariates=list(trt = 0, female = dat$female,
                                                          sofa = dat.boot$sofa, age = dat.boot$age))$`Expected percentages`[11, "State 1"]/100
            rd[i] <- p1 - p0

        }}
    rd
}

#### generate design matrices ##########################################################################################
#cross-sectional models
design_matrix_glm   <- model.matrix(mod.glm)
coefs_glm           <- list(coefs = coef(mod.glm), vcovs = vcov(mod.glm))
design_matrix_vgam1 <- model.matrix(mod.vgam1)
coefs_vgam1         <- list(coefs = coef(mod.vgam1), vcovs = vcov(mod.vgam1))
design_matrix_vgam2 <- model.matrix(mod.vgam2)
coefs_vgam2         <- list(coefs = coef(mod.vgam2), vcovs = vcov(mod.vgam2))
design_matrix_vgam3 <- model.matrix(mod.vgam3)
coefs_vgam3         <- list(coefs = coef(mod.vgam3), vcovs = vcov(mod.vgam3))

# longitudinal models
XMat        <- with(dat, cbind(trt, ns(logtime, 3),
                               trt*ns(logtime, 3), sofa, age, female))
XMat        <- XMat[1:nrow(XMat) %% 28 == 0,]
XMat.noX    <- with(dat, cbind(ns(logtime, 3),
                               trt*ns(logtime, 3), sofa, age, female))
XMat.noX    <- XMat.noX[1:nrow(XMat.noX) %% 28 == 0,]

XMat.ppo2    <- with(dat, cbind(ns(logtime, 3), ns(logtime, 3), trt*ns(logtime, 3)))
XMat.ppo2    <- XMat.ppo2[1:nrow(XMat.ppo2) %% 28 == 0,]

# bmtm
design_matrix_bmtm <- cbind(1, XMat)
colnames(design_matrix_bmtm ) <- c("Intercept", "trt", "t1", "t2", "t3", "trt*t1", "trt*t2", "trt*t3",
                                   "sofa", "age", "female")
coefs_bmtm         <- list(coefs = c(mtlv$alpha, mtlv$beta), vcovs = mtlv$vcov)
# bmtm.noX
design_matrix_bmtm.noX <- cbind(1, XMat.noX)
colnames(design_matrix_bmtm.noX) <- c("Intercept", "t1", "t2", "t3", "trt*t1", "trt*t2", "trt*t3",
                                      "sofa", "age", "female")
coefs_bmtm.noX         <- list(coefs = c(mtlv.noX$alpha, mtlv.noX$beta), vcovs = mtlv.noX$vcov)

# omtm1 and omtm3: the design matrix for the mortality part is the same as the bmtm
design_matrix_omtm1     <- design_matrix_bmtm
coefs_omtm1             <- list(coefs = c(omtm1$alpha[1], omtm1$beta[1:10]), vcovs = omtm1$vcov[c(1, 4:13),c(1, 4:13)])
design_matrix_omtm1.noX <- design_matrix_bmtm.noX
coefs_omtm1.noX         <- list(coefs = c(omtm1.noX$alpha[1], omtm1.noX$beta[1:9]),
                                vcovs = omtm1.noX$vcov[c(1, 4:12),c(1, 4:12)])
design_matrix_omtm3     <- design_matrix_bmtm
coefs_omtm3             <- list(coefs = c(omtm3$alpha[1], omtm3$beta[1:10]),
                                vcovs = omtm3$vcov[c(1, 4:13),c(1, 4:13)])
design_matrix_omtm3.noX <- design_matrix_bmtm.noX
coefs_omtm3.noX         <- list(coefs = c(omtm3.noX$alpha[1], omtm3.noX$beta[1:9]),
                                vcovs = omtm3.noX$vcov[c(1, 4:12),c(1, 4:12)])

design_matrix_omtm2     <- cbind(design_matrix_bmtm, XMat.ppo2[,c(1:3, 7:9)])
coefs_omtm2             <- list(coefs = c(omtm2$alpha[3], omtm2$beta[c(1:10, 14:19)]),
                                vcovs = omtm2$vcov[c(3, 4:13, 17:22),c(3, 4:13, 17:22)])
colnames(design_matrix_omtm2 ) <- c("Intercept", "trt", "t1", "t2", "t3", "trt*t1", "trt*t2", "trt*t3",
                                    "sofa", "age", "female",
                                    "t1:1", "t1:2", "t1:3",
                                    "trt*t1:1", "trt*t2:1", "trt*t3:1")
design_matrix_omtm2.noX <- cbind(design_matrix_bmtm.noX, XMat.ppo2[,c(1:3, 7:9)])
coefs_omtm2.noX             <- list(coefs = c(omtm2.noX$alpha[3], omtm2.noX$beta[c(1:9, 13:18)]),
                                    vcovs = omtm2.noX$vcov[c(3, 4:12, 16:21),c(3, 4:12, 16:21)])
colnames(design_matrix_omtm2.noX) <- c("Intercept", "t1", "t2", "t3", "trt*t1", "trt*t2", "trt*t3",
                                       "sofa", "age", "female",
                                       "t1:1", "t1:2", "t1:3",
                                       "trt*t1:1", "trt*t2:1", "trt*t3:1")

# survival
design_matrix_surv <- model.matrix(survs)
coefs_surv         <- list(coefs = survs$coefficients, vcovs = vcov(survs))







rd.glm <- rd.vgam1 <- rd.vgam2 <- rd.vgam3 <- rd.bmtm <- rd.bmtm.noX <- rd.omtm1 <-
    rd.omtm2 <- rd.omtm2.noX <- rd.omtm1.noX <- rd.omtm3 <- rd.omtm3.noX <-
    rd.survival <- c()


set.seed(91)

for(i in 1:5000){
    print(i)

    rd.glm[i]        <- rd_computation(coefficient = coefs_glm,  design_matrix = design_matrix_glm,
                                       method = "glm")
    rd.vgam1[i]      <- rd_computation(coefficient = coefs_vgam1,  design_matrix = design_matrix_vgam1,
                                       method = "vgam1")
    rd.vgam2[i]      <- rd_computation(coefficient = coefs_vgam2,  design_matrix = design_matrix_vgam2,
                                       method = "vgam2")
    rd.vgam3[i]      <- rd_computation(coefficient = coefs_vgam3,  design_matrix = design_matrix_vgam3,
                                       method = "vgam3")
    rd.bmtm[i]       <- rd_computation(coefficient = coefs_bmtm,  design_matrix = design_matrix_bmtm,
                                       method = "bmtm")
    rd.bmtm.noX[i]   <- rd_computation(coefficient = coefs_bmtm.noX,  design_matrix = design_matrix_bmtm.noX,
                                       method = "bmtm.noX")
    rd.omtm1[i]      <- rd_computation(coefficient = coefs_omtm1,  design_matrix = design_matrix_omtm1,
                                       method = "omtm1")
    rd.omtm1.noX[i]  <- rd_computation(coefficient = coefs_omtm1.noX,  design_matrix = design_matrix_omtm1.noX,
                                       method = "omtm1.noX")
    rd.omtm2[i]      <- -rd_computation(coefficient = coefs_omtm2,  design_matrix = design_matrix_omtm2,
                                        method = "omtm2")
    rd.omtm2.noX[i]  <- -rd_computation(coefficient = coefs_omtm2.noX,  design_matrix = design_matrix_omtm2.noX,
                                        method = "omtm2.noX")
    rd.omtm3[i]      <- rd_computation(coefficient = coefs_omtm3,  design_matrix = design_matrix_omtm3,
                                       method = "omtm3")
    rd.omtm3.noX[i]  <- rd_computation(coefficient = coefs_omtm3.noX,  design_matrix = design_matrix_omtm3.noX,
                                       method = "omtm3.noX")
    rd.survival[i]  <- rd_computation(coefficient = coefs_surv,  design_matrix = design_matrix_surv,
                                      method = "survival")
}

dat.final <- data.frame(Method = c("GLMB", "BMTM", "BMTM-C",
                                   "GLMO1", "GLMO2","GLMO3",
                                   "OMTM1", "OMTM1-C",
                                   "OMTM2", "OMTM2-C",
                                   "OMTM3", "OMTM3-C", "COXPH"),
                        estimated_rd = c(mean(rd.glm), mean(rd.bmtm), mean(rd.bmtm.noX),
                                         mean(rd.vgam1),mean(rd.vgam2), mean(rd.vgam3),
                                         mean(rd.omtm1), mean(rd.omtm1.noX),
                                         mean(rd.omtm2), mean(rd.omtm2.noX),
                                         mean(rd.omtm3), mean(rd.omtm3.noX), mean(rd.survival)
                        ),

                        se_rd = c(sd(rd.glm),   sd(rd.bmtm), sd(rd.bmtm.noX),
                                  sd(rd.vgam1), sd(rd.vgam2), sd(rd.vgam3),
                                  sd(rd.omtm1), sd(rd.omtm1.noX),
                                  sd(rd.omtm2), sd(rd.omtm2.noX),
                                  sd(rd.omtm3), sd(rd.omtm3.noX), sd(rd.survival)
                        ))
dat.final


# multistate model
mstate <- boot_rd(method = "multistate", nboot = 5000)
dat.final <- rbind(dat.final, data.frame(Method = "MULTISTATE", estimated_rd = mean(mstate), se_rd = sd(mstate)))
write.csv(dat.final , "./BinaryOrdinal28days/risk_difference.csv", row.names = FALSE)

