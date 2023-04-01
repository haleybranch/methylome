####### 
library(tidyverse)
library(lmerTest)
library(lme4)
library(lmtest) # for likelihood ratio tests

dat <- read.csv("S2_ft.csv")
s0 <- read.csv("s2_s0_ids.csv")

### merge dataframes together
dat2 <- inner_join(dat,s0)


## merge S0 and S1 treatments
dat2 <- dat2 %>%
  rowwise()%>%
  mutate(TRTHIS=paste(S0_Treatment,S1_treatment,sep=""))

dat2$FT <- as.numeric(dat2$FT)
dat2_filtered <- dat2 %>%
  filter(FT>0 & FT<1000 & !is.na(FT))

mod1 <- lmer(FT ~ Region * PrePeak * TRTHIS * Treatment + 
  (1|Site) + (1|S1_ID) + (1|S0_ID) + (1|REP), control=lmerControl(optimizer = "bobyqa", optCtrl=list(maxfun=100000)), data=dat2_filtered)

anova(mod1)

### drop region prepeak s1_treatment s0_treatment treatment

mod2 <- lmer(FT ~ Region*PrePeak*S1_treatment*S0_Treatment + 
               Region*PrePeak*S1_treatment*Treatment +
               Region*S1_treatment*S0_Treatment*Treatment +
               PrePeak*S1_treatment*S0_Treatment*Treatment +
               (1|Site) + (1|S1_ID) + (1|S0_ID) + (1|REP), control=lmerControl(optimizer = "bobyqa", optCtrl=list(maxfun=100000)), data=dat2_filtered)

lrtest(mod1, mod2) # accept 4-way model

########
datdouble <- rbind(dat2_filtered,dat2_filtered)

mod1 <- lmer(FT ~ Region * PrePeak * TRTHIS * Treatment + 
               (1|Site) + (1|S1_ID) + (1|S0_ID), control=lmerControl(optimizer = "bobyqa", optCtrl=list(maxfun=100000)), data=datdouble)
anova(mod1)




##### gt

library(tidyverse)
library(lmerTest)
library(lme4)
library(lmtest) # for likelihood ratio tests

dat <- read.csv("S2_gt.csv")
s0 <- read.csv("s2_s0_ids.csv")

### merge dataframes together
dat2 <- inner_join(dat,s0)


## merge S0 and S1 treatments
dat2 <- dat2 %>%
  rowwise()%>%
  mutate(TRTHIS=paste(S0_Treatment,S1_treatment,sep=""))

dat2$days_to_germ <- as.numeric(dat2$days_to_germ)
dat2_filtered <- na.omit(dat2)


mod1 <- lmer(days_to_germ ~ Region * PrePeak * TRTHIS + 
               (1|Site) + (1|S1_ID) + (1|S0_ID) + (1|Rep), control=lmerControl(optimizer = "bobyqa", optCtrl=list(maxfun=100000)), data=dat2_filtered)

anova(mod1)

# drop 3way
mod2 <- lmer(days_to_germ ~ Region*PrePeak + TRTHIS*Region +
                       PrePeak*TRTHIS + (1|S1_ID) + (1|S0_ID) + (1|Rep) + (1|Site), 
                     control=lmerControl(optimizer = "bobyqa", optCtrl=list(maxfun=100000)), data=dat2_filtered)

lrtest(mod1, mod2) # accept 3-way model


library(ggplot2)
library(visreg)

vis_south<-visreg(mod1, xvar="PrePeak", by="TRTHIS", cond=list(Treatment="South"))

vis_north<-visreg(mod1, xvar="PrePeak", by="TRTHIS", cond=list(Treatment="North"))

Res_south<-vis_south$res ; Res_north<-vis_north$res # Extract residuals
Res_all<-rbind(Res_south, Res_north) #Row bind wet and dry residuals into one data frame


