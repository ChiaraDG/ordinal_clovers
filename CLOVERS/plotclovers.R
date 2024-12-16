############################################################################################################################
############## Plot CLOVERS Data ###########################################################################################
############################################################################################################################

###### load the data #######################################################################################################
library(tidyverse)
library(viridis)
library(OMTM1)
library(knitr)
library(VGAM)
library(forestploter)
library(cowplot)

# load the outcome
dat           <- read.csv("/Users/chiaradigravio/Documents/PhD/Collaborations/CLOVERS/Data/dailystateICU28daysmortality_oneabssorbing.csv")
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

#### Transition probabilities
mat <- Calc.TransProbs(dat$DailyState, id = dat$ReferenceID, digs = 4)$Tx.Mtx.Col2Row
kable(mat, digits = 3)

###### Figure 1 #######################################################################################################

# panel A: histogram
p1 <- dat %>%
  dplyr::group_by(Day, DailyState) %>% dplyr::summarize(num = sum(!is.na(DailyState)), .groups = "drop") %>%
  na.omit() %>% group_by(Day) %>% mutate(tot = sum(num), prop = num/tot) %>%
  ggplot(aes(x = Day, y = prop, fill = DailyState)) +
  geom_bar(position="stack", stat="identity") +  scale_fill_viridis(discrete = T) +
  theme_minimal() +  theme(legend.position = "bottom") +
  labs(x = "Days from Randomization", y = "Proportion of Patients",
       title = "Proportion of Subjects in Each State\nin the First 28 Days Post-Randomization") +
  scale_x_continuous(breaks = seq(0, 28, by = 2))


# panel B: cumulative log-odds
lo <- function(x) log(mean(x)/(1-mean(x)))

plot.states.lo <- function(data, Yvar, timevar){

    dat1       = data
    dat1$t     = dat1[,timevar]
    lid        = length(unique(dat1$id))

    l1 = tapply(dat1[,Yvar]<=1, dat1$t, lo)
    l2 = tapply(dat1[,Yvar]<=2, dat1$t, lo)
    l3 = tapply(dat1[,Yvar]<=3, dat1$t, lo)
    times = sort(unique(dat1$t))

    out = data.frame(times = times, l1 = l1, l2 = l2, l3 = l3)

    out = pivot_longer(out, cols = l1:l3)

    out
}

tmp <- plot.states.lo(data = subset(dat, Day > 0), Yvar = "DailyState2", timevar = "Day")
col.vir=viridis(20)

p2 <- ggplot(tmp, aes(x = times, y = value, group = name)) +
    geom_ribbon(data = subset(tmp, name == "l3"), aes(ymin = -Inf, ymax = value), fill = col.vir[8]) +
    geom_ribbon(data = subset(tmp, name == "l2"), aes(ymin = -Inf, ymax = value), fill = col.vir[15]) +
    geom_ribbon(data = subset(tmp, name == "l3"), aes(ymax = 4.7,  ymin = value), fill = col.vir[1]) +
    geom_ribbon(data = subset(tmp, name == "l1"), aes(ymin = -Inf, ymax = value), fill = col.vir[20]) +
    geom_line(linewidth = 0.7) + theme_minimal() +
    scale_x_continuous(breaks = seq(0, 28, by = 2)) +
    labs(x = "Days from Randomization", y = "Cumulative Log-Odds",
         title = "Empirical Cumulative Log-Odds\nin the First 28 Days Post-Randomization")

legend_b <- get_legend(p1 + theme(legend.position="bottom"))

prow <- plot_grid( p1 + theme(legend.position="none"),
                   p2,
                   align = 'vh',
                   labels = c("A", "B"),
                   hjust = -1,
                   nrow = 1
)

ptot <- plot_grid(prow, legend_b, ncol = 1, rel_heights = c(1.2, .1))
ptot

###### Figure 3 #######################################################################################################

# transform an ordinal data to a binary data
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
omtm3    <- readRDS(paste0(dir, "omtm3.RData"))
omtm3.noX <- readRDS(paste0(dir, "omtm3_noX.RData"))



basis <- ns(dat$logtime, 3)
prdns <- predict(basis, 1)

est.coef <- c(
  # GLM
  coef(mod.glm)[2],
  # BMTM
  cbind(1, 1*prdns)%*%(mtlv$beta[c(1, 5:7)]),
  # BMTM-C
  cbind(1*prdns)%*%(mtlv.noX$beta[c(4:6)]),
  # GLMO1, GLMO2, GLMO3
  coef(mod.vgam1)[4], coef(mod.vgam2)[4], coef(mod.vgam3)[4],
  # OMTM1
  cbind(1, 1*prdns)%*%(omtm1$beta[c(1, 5:7)]),
  # OMTM proportional odds-C
  cbind(1*prdns)%*%(omtm1.noX$beta[c(4:6)]),
  # OMTM2
  -cbind(1, 1*prdns, 1*prdns)%*%(omtm2$beta[c(1, 5:7, 17:19)]),
  # OMTM2-C
  -cbind(1*prdns, 1*prdns)%*%(omtm2.noX$beta[c(4:6, 16:18)]),
  # OMTM3
  cbind(1, 1*prdns)%*%(omtm3$beta[c(1, 5:7)]),
  # OMTM3-C
  cbind(1*prdns)%*%(omtm3.noX$beta[c(4:6)])
)

mat.ppo.deathonly <- cbind(
  omtm3$vcov[c(4, 8:10, 20:22), c(4, 8:10, 20:22)])

mat.ppo.deathonly.noX <- cbind(
  omtm3.noX$vcov[c(7:9, 19:21), c(7:9, 19:21)])

#mat.ppo.deathonly[5:7, 5:7] <- -mat.ppo.deathonly[5:7, 5:7]
se.coef  <- sqrt(c(vcov(mod.glm)[2,2],
                   # BMTM
                   cbind(1, 1*prdns)%*%mtlv$vcov[c(2, 6:8),c(2, 6:8)]%*%t(cbind(1, 1*prdns)),
                   cbind(1*prdns)%*%mtlv.noX$vcov[c(5:7),c(5:7)]%*%t(cbind(1*prdns)),
                   # VGAM
                   vcov(mod.vgam1)[4,4], vcov(mod.vgam2)[4,4], vcov(mod.vgam3)[4,4],
                   # OMTM1
                   cbind(1, 1*prdns)%*%omtm1$vcov[c(4, 8:10),c(4, 8:10)]%*%t(cbind(1, 1*prdns)),
                   cbind(1*prdns)%*%omtm1.noX$vcov[c(7:9),c(7:9)]%*%t(cbind(1*prdns)),
                   # OMTM2
                   cbind(1, 1*prdns, 1*prdns)%*%mat.ppo.deathonly%*%t(cbind(1, 1*prdns, 1*prdns)),
                   cbind(1*prdns, 1*prdns)%*%mat.ppo.deathonly.noX%*%t(cbind(1*prdns, 1*prdns)),
                   # OMTM3
                   cbind(1, 1*prdns)%*%omtm3$vcov[c(4, 8:10),c(4, 8:10)]%*%t(cbind(1, 1*prdns)),
                   cbind(1*prdns)%*%omtm3.noX$vcov[c(7:9),c(7:9)]%*%t(cbind(1*prdns))))

res <- data.frame(Method = c("GLMB", "BMTM",
                             "BMTM-C",
                             "GLMO1", "GLMO2",
                             "GLMO3",
                             "OMTM1", "OMTM1-C",
                             "OMTM2", "OMTM2-C",
                             "OMTM3", "OMTM3-C"),
                  coefs = est.coef, ses = se.coef)

res$lb <- res$coefs - qnorm(0.975)*res$ses
res$ub <- res$coefs + qnorm(0.975)*res$ses

res$`log(OR) (95% CI)` <- ifelse(is.na(res$ses), "",
                                 sprintf("%.3f (%.3f, %.3f)",
                                         res$coefs, res$lb, res$ub))
res$` ` <- paste(rep(" ", 20), collapse = " ")
res$`CI Width` <- sprintf("%.3f", res$ub - res$lb)

tm <- forest_theme(base_size = 12,
                   # Reference line width/type/color
                   refline_lwd = 1,
                   refline_lty = "dashed",
                   refline_col = "grey20",
                   core = list(padding = unit(c(3, 1), "cm")))

# risk difference
rd_dat <- read.csv("~/Documents/PhD/Collaborations/CLOVERS/BinaryOrdinal28days/risk_difference.csv")

rd_dat$Method[rd_dat$Method == "MULTISTATE"] <- "MSM"


rd_dat$lb <- rd_dat$estimated_rd - qnorm(0.975)*rd_dat$se_rd
rd_dat$ub <- rd_dat$estimated_rd + qnorm(0.975)*rd_dat$se_rd

rd_dat$`RD (95% CI)` <- ifelse(is.na(rd_dat$se_rd), "",
                               sprintf("%.3f (%.3f, %.3f)",
                                       rd_dat$estimated_rd, rd_dat$lb, rd_dat$ub))
rd_dat$` ` <- paste(rep(" ", 20), collapse = " ")
rd_dat$`CI Width` <- sprintf("%.3f", rd_dat$ub - rd_dat$lb)

# put odds ratio and risk difference in tha same plot
tmp <- full_join(res, rd_dat, by = "Method")
names(tmp) <- c("Method", "log_or", "ses_or", "lb_or", "ub_or",
                "log(OR) (95% CI)", "", "log(OR)\nCI Width",
                "rd", "ses_rd", "lb_rd", "ub_rd", "RD (95% CI)", " ", "RD\nCI Width")

p <- forest(tmp[,c(1, 6:8, 13:15)],
            est = list(tmp$log_or,
                       tmp$rd),
            lower = list(tmp$lb_or,
                         tmp$lb_rd),
            upper = list(tmp$ub_or,
                         tmp$ub_rd),
            ci_column = c(3, 6),
            ref_line = c(0, 0), theme = tm,
            xlab = c("Log-Odds Ratio", "Risk Difference"))
p

ggsave("CloversResultsAll.png", p,  units = "cm", width = 50, height = 25)
