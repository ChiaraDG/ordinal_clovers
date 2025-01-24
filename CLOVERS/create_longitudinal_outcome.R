library(tidyverse)
library(viridis)

dat <- read.csv("Data/dataset.csv") # baseline information + days from randomization to an event

mvent <- unique(dat[,c("ReferenceID", "venthx_intdt", "venthx_ventyn", "venthx_studyint",
                       "venthx_randint", "venthx_startdt", "venthx_lastdt", "ardshx_randyn", "ardshx_randd7yn",
                       names(dat)[grep("^icu", names(dat))])])
id.todrop <- mvent %>% filter(mvent$icuhx_icuyn == "") %>%
    select(ReferenceID) %>% pull()

dat.small <- unique(dat[c("ReferenceID", "venthx_intdt", "venthx_randint", "venthx_startdt", "ardshx_randd7yn",
                          "venthx_lastdt", "term_dischdt", "deathdt",
                          "term_wddt", names(dat)[grep("^icu", names(dat))])])

id.todrop2 <- dat.small %>% filter(ardshx_randd7yn == "Yes" & !is.na(deathdt)) %>% select(ReferenceID) %>% pull()
todrop <- unique(c(id.todrop, id.todrop2))
# 34 people total should be dropped
dat.small <- dat.small[!c(dat.small$ReferenceID %in% todrop), ]

############ CREATE DAILY STATE DATA ############
# If someone died or they ware discharged after day 28 set them as NA
dat.small$term_dischdt <- ifelse(dat.small$term_dischdt > 28, NA, dat.small$term_dischdt)
dat.small$deathdt      <- ifelse(dat.small$deathdt > 28, NA, dat.small$deathdt)
# If someone has enter from ICU greater than 28 set it to NA
dat.small$icuhx_admitdt_r1  <- ifelse(dat.small$icuhx_admitdt_r1 > 28, NA, dat.small$icuhx_admitdt_r1)
dat.small$icuhx_admitdt_r2  <- ifelse(dat.small$icuhx_admitdt_r2 > 28, NA, dat.small$icuhx_admitdt_r2)
dat.small$icuhx_admitdt_r3  <- ifelse(dat.small$icuhx_admitdt_r3 > 28, NA, dat.small$icuhx_admitdt_r3)
dat.small$icuhx_admitdt_r4  <- ifelse(dat.small$icuhx_admitdt_r4 > 28, NA, dat.small$icuhx_admitdt_r4)
dat.small$icuhx_admitdt_r5  <- ifelse(dat.small$icuhx_admitdt_r5 > 28, NA, dat.small$icuhx_admitdt_r5)
# If someone has exits from ICU greater than 28 set it to 28
dat.small$icuhx_dischargedt_r1  <- ifelse(dat.small$icuhx_dischargedt_r1 > 28, 28, dat.small$icuhx_dischargedt_r1)
dat.small$icuhx_dischargedt_r2  <- ifelse(dat.small$icuhx_dischargedt_r2 > 28, 28, dat.small$icuhx_dischargedt_r2)
dat.small$icuhx_dischargedt_r3  <- ifelse(dat.small$icuhx_dischargedt_r3 > 28, 28, dat.small$icuhx_dischargedt_r3)
dat.small$icuhx_dischargedt_r4  <- ifelse(dat.small$icuhx_dischargedt_r4 > 28, 28, dat.small$icuhx_dischargedt_r4)
dat.small$icuhx_dischargedt_r5  <- ifelse(dat.small$icuhx_dischargedt_r5 > 28, 28, dat.small$icuhx_dischargedt_r5)
# If someone has an enter date in the ICU but not the exit, assume they never exit the ICU
# Example ID: D01-01333
dat.small$icuhx_dischargedt_r1  <- ifelse(!is.na(dat.small$icuhx_admitdt_r1) & is.na(dat.small$icuhx_dischargedt_r1), 28,
                                          dat.small$icuhx_dischargedt_r1)
dat.small$icuhx_dischargedt_r2  <- ifelse(!is.na(dat.small$icuhx_admitdt_r2) & is.na(dat.small$icuhx_dischargedt_r2), 28,
                                          dat.small$icuhx_dischargedt_r2)
dat.small$icuhx_dischargedt_r3  <- ifelse(!is.na(dat.small$icuhx_admitdt_r3) & is.na(dat.small$icuhx_dischargedt_r3), 28,
                                          dat.small$icuhx_dischargedt_r3)
dat.small$icuhx_dischargedt_r4  <- ifelse(!is.na(dat.small$icuhx_admitdt_r4) & is.na(dat.small$icuhx_dischargedt_r4), 28,
                                          dat.small$icuhx_dischargedt_r4)
dat.small$icuhx_dischargedt_r5  <- ifelse(!is.na(dat.small$icuhx_admitdt_r5) & is.na(dat.small$icuhx_dischargedt_r5), 28,
                                          dat.small$icuhx_dischargedt_r5)


head(dat.small)
# people who died after discharge
dim(dat.small[!is.na(dat.small$term_dischdt) & !is.na(dat.small$deathdt),])
# 68 people have a discharge date later than a death date

# create a dataset where each person is observed 29 times.
# set state = to hospital for each subject
dat.long <- dat.small[rep(seq_len(1529), each = 29), ]
dat.long$Day <- 0:28

create.state       <- function(ID, d.long = dat.long){
    tmp.28   = subset(d.long,  ReferenceID == ID)
    # set everyone at home to start with
    DayState = rep("Hospital", 29)
    # set those who died
    for(i in 1:29){
        DayState[i] <- ifelse(!is.na(tmp.28$deathdt[1]) & tmp.28$deathdt[1] <= tmp.28$Day[i],
                              "Death", "Hospital")
    }
    # set up those who are in the ICU
    # Condition 1: at some point patients where in the ICU
    condition1 <- (tmp.28$icuhx_icuyn[1] == "Yes" & !is.na(tmp.28$icuhx_icuyn[1]))
    # condition 2: Patients who were in the ICU at least one time
    condition2 <- (condition1 & !is.na(tmp.28$icuhx_admitdt_r1[1]))
    # condition 3: Patients who were in the ICU at least two times
    condition3 <- (condition1 & !is.na(tmp.28$icuhx_admitdt_r2[1]))
    # condition 4: Patients who were in the ICU at least three time
    condition4 <- (condition1 & !is.na(tmp.28$icuhx_admitdt_r3[1]))
    # condition 5: Patients who were in the ICU at least four times
    condition5 <- (condition1 & !is.na(tmp.28$icuhx_admitdt_r4[1]))
    # condition 6: Patients who were in the ICU five times
    condition6 <- (condition1 & !is.na(tmp.28$icuhx_admitdt_r5[1]))

    if(condition1 & condition2){
        for(i in 1:29){
            DayState[i] <- ifelse(tmp.28$icuhx_admitdt_r1[1] <= tmp.28$Day[i] &
                                      tmp.28$icuhx_dischargedt_r1[1] >= tmp.28$Day[i] &
                                      DayState[i] != "Death",
                                  "ICU", DayState[i])
        }
    }

    if(condition1 & condition3){
        for(i in 1:29){
            DayState[i] <- ifelse(tmp.28$icuhx_admitdt_r2[1] <= tmp.28$Day[i] &
                                      tmp.28$icuhx_dischargedt_r2[1] >= tmp.28$Day[i] &
                                      DayState[i] != "Death",
                                  "ICU", DayState[i])
        }
    }

    if(condition1 & condition4){
        for(i in 1:29){
            DayState[i] <- ifelse(tmp.28$icuhx_admitdt_r3[1] <= tmp.28$Day[i] &
                                      tmp.28$icuhx_dischargedt_r3[1] >= tmp.28$Day[i] &
                                      DayState[i] != "Death",
                                  "ICU", DayState[i])
        }
    }

    if(condition1 & condition5){
        for(i in 1:29){
            DayState[i] <- ifelse((tmp.28$icuhx_admitdt_r4[1] <= tmp.28$Day[i]) &
                                      tmp.28$icuhx_dischargedt_r4[1] >= tmp.28$Day[i] &
                                      DayState[i] != "Death",
                                  "ICU", DayState[i])
        }
    }

    if(condition1 & condition6){
        for(i in 1:29){
            DayState[i] <- ifelse(tmp.28$icuhx_admitdt_r5[1] <= tmp.28$Day[i] &
                                      tmp.28$icuhx_dischargedt_r5[1] >= tmp.28$Day[i] &
                                      DayState[i] != "Death",
                                  "ICU", DayState[i])
        }
    }

    # set up those who were discharged. If someone is discharged the same day
    # they were taken off ICU those count as ICU days
    condition7 <- !is.na(tmp.28$term_dischdt[1])
    if(condition7){
        for(i in 1:29){
            DayState[i] <- ifelse((DayState[i] != "Death" & DayState[i] != "ICU") &
                                      tmp.28$term_dischdt[i] < tmp.28$Day[i], "Discharged",
                                  DayState[i])
        }
    }

    return(DayState)
}

tmp                 <- lapply(X = unique(dat.long$ReferenceID), FUN = create.state)
dat.long$DailyState <- unlist(tmp)
dat.long$DailyState <- ordered(dat.long$DailyState,
                               levels = c("Discharged", "Hospital", "ICU",
                                          "Death"))
dat.long            <- subset(dat.long, select = c("ReferenceID", "Day", "DailyState"))
initial.state       <- subset(dat.long, Day == 0, select = c("ReferenceID", "DailyState"))
names(initial.state) <- c("ReferenceID", "InitialState")
dat.long            <- inner_join(dat.long, initial.state, by = "ReferenceID")
write.csv(dat.long, "./Data/dailystateICU28daysmortality_oneabssorbing.csv", row.names = FALSE)

