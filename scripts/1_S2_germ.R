############## Germination date #################
# Written by Haley Branch August 18, 2022

####### load in libraries
library(tidyverse)
library(lmerTest)
library(lme4)
library(ggplot2)
library(visreg)
library(emmeans)
library(multtest)
library(effectsize)


# read in data
dat_pre <- read.csv("Data/ft_gt_pre.csv")
dat_peak <- read.csv("Data/ft_gt_peak.csv")

#filter out one of the treatments, because otherwise the data is doubled
#germination data does not have a current treatment

dat_pre <- filter(dat_pre, Treatment == "Wet")
dat_peak <- filter(dat_peak, Treatment == "Wet")
dat <- rbind(dat_pre, dat_peak)

#convert variables
dat$ID<- as.character(dat$ID)
dat$FT<- as.numeric(dat$FT)
dat$FD<- as.numeric(dat$FD)
dat$days_to_germ<- as.numeric(dat$days_to_germ)

# break up the data by region and treatment 
dat_north <- filter(dat, Region == "North")
dat_south <- filter(dat, Region == "South")

dat_north <- dat_north %>% filter(!is.na(days_to_germ))
dat_south <- dat_south %>% filter(!is.na(days_to_germ))

#grab trait means
means <- dat_north %>%
  group_by(TRTHIS) %>%
  summarise_at(vars(days_to_germ), list(name = mean))
means

#set up model pre
mod_north <- lmer(days_to_germ ~ Rep + Site + PrePeak*TRTHIS + (1|S0_ID), data=dat_north)

lattice::qqmath(mod_north) # qqplot

#log transform the data
dat_north$germ_log <- log(dat_north$days_to_germ)

mod_north_log <- lmer(germ_log ~ Rep + Site + PrePeak*TRTHIS + (1|S0_ID), data=dat_north)

lattice::qqmath(mod_north_log) # qqplot


hist(resid(mod_north_log)) #histogram


plot(mod_north_log, rstudent(.) ~ hatvalues(.)) #residuals vs leverage

#anova
anova_north <- anova(mod_north_log) 
anova_north
#trthist p = 0.0003 
#write.csv(anova_north, "Results/germ_north_anova.csv")


#grab effect sizes of the fixed factors
#es = difference between groups/standard deviation of the variance 
effsize_north <- omega_squared(mod_north_log)
effsize_north
#TRTHIS = 0.04
#write.csv(effsize_north, "Results/germ_north_effectsize.csv")


### create planned contrasts

emm1 = emmeans(mod_north_log, specs = ~ PrePeak*TRTHIS)
emm1


peak_DD = c(1, 0, 0, 0, 0, 0, 0, 0)
pre_DD = c(0, 1, 0, 0, 0, 0, 0, 0)
peak_DW = c(0, 0, 1, 0, 0, 0, 0, 0)
pre_DW = c(0, 0, 0, 1, 0, 0, 0, 0)
peak_WD = c(0, 0, 0, 0, 1, 0, 0, 0)
pre_WD = c(0, 0, 0, 0, 0, 1, 0, 0)
peak_WW = c(0, 0, 0, 0, 0, 0, 1, 0)
pre_WW = c(0, 0, 0, 0, 0, 0, 0, 1)



contrasts <- list("pre_DD-DW" = pre_DD - pre_DW, "pre_DD-WD" = pre_DD - pre_WD,
                  "pre_DD-WW" = pre_DD - pre_WW, 
                  "pre_DW-WW" = pre_DW - pre_WW, "pre_WD-WW" = pre_WD - pre_WW,
                  "peak_DD-DW" = peak_DD - peak_DW, "peak_DD-WD" = peak_DD - peak_WD,
                  "peak_DD-WW" = peak_DD - peak_WW, 
                  "peak_DW-WW" = peak_DW - peak_WW, "peak_WD-WW" = peak_WD - peak_WW,
                  
                  "pre-peak_DD" = pre_DD - peak_DD, "pre-peak_DW" = pre_DW - peak_DW,
                  "pre-peak_WD" = pre_WD - peak_WD, "pre-peak_WW" = pre_WW - peak_WW)


#grab the raw pvalues without a posthoc adjustment
contrast.amax <- contrast(emm1, method = contrasts, adjust="none") 
summ.contrast <- summary(contrast.amax)

#use a pvalue correction for the multiple tests 
adjust.p <- mt.rawp2adjp(summ.contrast$p.value, proc=c("BH"), alpha = 0.05, na.rm = FALSE)
#adjust.p <- as.data.frame(adjust.p$adjp)
#adjust.p <- filter(adjust.p, rawp > 0.0001)
#realign the adjusted pvalues with the contrast names
adjust.p<-as.data.frame(adjust.p$adjp)
adjust.p<-distinct(adjust.p)
adjust.p.df<- adjust.p %>%
  rename(p.value=rawp)%>%
  inner_join(summ.contrast)
adjust.p.df
#write.csv(adjust.p.df, "Results/germ_north_planned.csv")
#pre DD-DW = 0.07 --> DD earlier than DW 
#peak DD-DW = 0.07 --> DD earlier than DW 


vis_germ_north<-visreg(mod_north_log, xvar="TRTHIS", by="PrePeak") 
Res_germ_north<-vis_germ_north$res # Extract residuals

Res_germ_north$PrePeak <- factor(Res_germ_north$PrePeak, level = c('Pre', 'Peak'))


plot_north <- ggplot(Res_germ_north, aes(TRTHIS, y=exp(visregRes)))+
  geom_violin(aes(), fill ="white", position = position_dodge(width = 1))+
  stat_summary(aes(), fun=mean, colour="black",position = position_dodge(width = 1))+
  stat_summary(aes(),fun.data=mean_se, geom="errorbar", colour="black", width=0.5,position = position_dodge(width = 1))+
  xlab("Treatment History") +
 facet_wrap(.~ PrePeak) +
  scale_y_continuous(name="Germination Time")+
  theme_classic()
plot_north <-plot_north + theme(
  axis.text.x = element_text(size=12, face="bold", angle=0,hjust=0.5),
  axis.text.y = element_text(size=15,face="bold"),
  axis.title.x = element_text(color="black", size=20, vjust = 0.5, face="bold"),
  axis.title.y = element_text(color="black", size=20,vjust = 2, face="bold",hjust=0.5))
plot_north <-plot_north + facet_wrap(.~ PrePeak) +
  theme(legend.title = element_blank(),legend.text = element_text(size=12,face="bold"),
        strip.background = element_blank(), strip.text.x=element_text(size=14,face="bold",hjust=0.05,vjust=-1.2))
plot_north
#save 6x4

################### SOUTH ############################################

mod_south <- lmer(days_to_germ ~ Rep + Site + PrePeak*TRTHIS + (1|S0_ID), data=dat_south)

lattice::qqmath(mod_south) # qqplot

#log transform the data
dat_south$germ_log <- log(dat_south$days_to_germ)

mod_south_log <- lmer(germ_log ~ Rep + Site + PrePeak*TRTHIS + (1|S0_ID), data=dat_south)

lattice::qqmath(mod_south_log) # qqplot


hist(resid(mod_south)) #histogram
hist(resid(mod_south_log))

plot(mod_south_log, rstudent(.) ~ hatvalues(.)) #residuals vs leverage

#anova
anova_south <- anova(mod_south_log) 
anova_south
#trthist p = 0.08
#year p = 0.06
#write.csv(anova_south, "Results/germsouth_anova.csv")


#grab effect sizes of the fixed factors
#es = difference between groups/standard deviation of the variance 
effsize_south <- omega_squared(mod_south_log)
effsize_south
#PrePeak = 0.14
#TRTHIS < 0.01
#write.csv(effsize_south, "Results/germsouth_effectsize.csv")


### create planned contrasts

emm1 = emmeans(mod_south_log, specs = ~ PrePeak*TRTHIS)
emm1


peak_DD = c(1, 0, 0, 0, 0, 0, 0, 0)
pre_DD = c(0, 1, 0, 0, 0, 0, 0, 0)
peak_DW = c(0, 0, 1, 0, 0, 0, 0, 0)
pre_DW = c(0, 0, 0, 1, 0, 0, 0, 0)
peak_WD = c(0, 0, 0, 0, 1, 0, 0, 0)
pre_WD = c(0, 0, 0, 0, 0, 1, 0, 0)
peak_WW = c(0, 0, 0, 0, 0, 0, 1, 0)
pre_WW = c(0, 0, 0, 0, 0, 0, 0, 1)



contrasts <- list("pre_DD-DW" = pre_DD - pre_DW, "pre_DD-WD" = pre_DD - pre_WD,
                  "pre_DD-WW" = pre_DD - pre_WW, 
                  "pre_DW-WW" = pre_DW - pre_WW, "pre_WD-WW" = pre_WD - pre_WW,
                  "peak_DD-DW" = peak_DD - peak_DW, "peak_DD-WD" = peak_DD - peak_WD,
                  "peak_DD-WW" = peak_DD - peak_WW, 
                  "peak_DW-WW" = peak_DW - peak_WW, "peak_WD-WW" = peak_WD - peak_WW,
                  
                  "pre-peak_DD" = pre_DD - peak_DD, "pre-peak_DW" = pre_DW - peak_DW,
                  "pre-peak_WD" = pre_WD - peak_WD, "pre-peak_WW" = pre_WW - peak_WW)


#grab the raw pvalues without a posthoc adjustment
contrast.amax <- contrast(emm1, method = contrasts, adjust="none") 
summ.contrast <- summary(contrast.amax)

#use a pvalue correction for the multiple tests 
adjust.p <- mt.rawp2adjp(summ.contrast$p.value, proc=c("BH"), alpha = 0.05, na.rm = FALSE)
#adjust.p <- as.data.frame(adjust.p$adjp)
#adjust.p <- filter(adjust.p, rawp > 0.0001)
#realign the adjusted pvalues with the contrast names
adjust.p<-as.data.frame(adjust.p$adjp)
adjust.p<-distinct(adjust.p)
adjust.p.df<- adjust.p %>%
  rename(p.value=rawp)%>%
  inner_join(summ.contrast)
adjust.p.df
#write.csv(adjust.p.df, "Results/germ_south_planned.csv")
#nothing

vis_germ_south<-visreg(mod_south_log, xvar="TRTHIS", by="PrePeak") 
Res_germ_south<-vis_germ_south$res # Extract residuals

Res_germ_south$PrePeak <- factor(Res_germ_south$PrePeak, level = c('Pre', 'Peak'))


plot_south <- ggplot(Res_germ_south, aes(TRTHIS, y=exp(visregRes)))+
  geom_violin(aes(), fill ="white", position = position_dodge(width = 1))+
  stat_summary(aes(), fun=mean, colour="black",position = position_dodge(width = 1))+
  stat_summary(aes(),fun.data=mean_se, geom="errorbar", colour="black", width=0.5,position = position_dodge(width = 1))+
  xlab("Treatment History") +
 facet_wrap(.~ PrePeak) +
  ylab("Germination Time") +
 # ymax(15)
  theme_classic()
plot_south <-plot_south + theme(
  axis.text.x = element_text(size=12, face="bold", angle=0,hjust=0.5),
  axis.text.y = element_text(size=15,face="bold"),
  axis.title.x = element_text(color="black", size=20, vjust = 0.5, face="bold"),
  axis.title.y = element_text(color="black", size=20,vjust = 2, face="bold",hjust=0.5))
plot_south <-plot_south +  facet_wrap(.~ PrePeak) +
  theme(legend.title = element_blank(),legend.text = element_text(size=12,face="bold"),
        strip.background = element_blank(), strip.text.x=element_text(size=14,face="bold",hjust=0.05,vjust=-1.2))
plot_south
#save 8x6