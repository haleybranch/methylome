####### sla
library(tidyverse)
library(lmerTest)
library(lme4)
library(ggplot2)
library(visreg)
library(emmeans)
library(multtest)
library(effectsize)
#s0

dat1<- read.csv("Data/s2_SLA_rep1.csv")
dat1$sla <-as.numeric(dat1$sla)
dat1 <- filter(dat1, !is.na(sla))
dat1$Rep <- "1"
dat1$X <- NULL
dat1$X.1 <- NULL
dat1$X.2 <- NULL
dat1$X.3 <- NULL

dat2<-read.csv("Data/s2_SLA_rep2.csv")
dat2$sla <-as.numeric(dat2$sla)
dat2 <- filter(dat2, !is.na(sla))
dat2$Rep <- "2"
dat <- rbind(dat1,dat2)


#s2 <- read.csv("Data/S2_IDs.csv")
s0 <- read.csv("Data/s2_s0_ids.csv")
trthis <- read_csv("Data/Trthis_IDs.csv")

### merge dataframes together
data <- merge(s0,trthis, by="S0_ID")
dat <- merge(dat,data, by="ID")

dat$Treatment = if_else(dat$Treatment == "Dry", "D", "W")

## merge S0 and S1 treatments
dat <- dat %>%
  rowwise()%>%
  mutate(Treatment_history=paste(S0_trt,S1_trt,Treatment,sep=""))

dat$Rep <- as.character(dat$Rep)

## Question: Is transgeneration plasticity present?
# break up the data by region and treatment 
dat_north_pre_site15 <- filter(dat, Region == "North", PrePeak == "Pre", Site =="S15")
dat_north_peak_site15 <- filter(dat, Region == "North", PrePeak == "Peak", Site =="S15")
dat_north_pre_site36 <- filter(dat, Region == "North", PrePeak == "Pre", Site =="S36")
dat_north_peak_site36 <- filter(dat, Region == "North", PrePeak == "Peak", Site =="S36")
#dat_south_pre <- filter(dat, Region == "South", PrePeak == "Pre", Site =="S15")
#dat_south_peak <- filter(dat, Region == "South",PrePeak == "Peak", Site =="S15")

#remove sla numbers greater than 900, since these are outliers 

dat_north_pre_site15 <- filter(dat_north_pre_site15, sla < 900)
dat_north_peak_site15 <- filter(dat_north_pre_site15, sla < 900)
dat_north_pre_site36 <- filter(dat_north_pre_site36, sla < 900)
dat_north_peak_site36 <- filter(dat_north_pre_site36, sla < 900)

################### NORTH_PRE ############################################
#### WET ######
mod_nw_pre_36 <- lm(sla~Rep + S0_ID*Treatment_history, data=dat_north_pre_site36)

#check for model violations
plot(mod_nw, type=c("p","smooth"), col.line=1) #fitted vs residuals

plot(mod_nw,
     sqrt(abs(resid(.)))~fitted(.),
     type=c("p","smooth"), col.line=1) #scale-location plot

lattice::qqmath(mod_nw) # qqplot

hist(resid(mod_nw)) #histogram

plot(mod_nw, rstudent(.) ~ hatvalues(.)) #residuals vs leverage

summary(mod_nw)

#anova
anova_nw <- anova(mod_nw) 
anova_nw
#no differences 
#write.csv(anova_nw, "Results/sla_nw_anova.csv")


#grab effect sizes of the fixed factors
#es = difference between groups/standard deviation of the variance 
effsize_north <- omega_squared(mod_nw)
effsize_north
#write.csv(effsize_north, "Results/sla_nw_effectsize.csv")

### create planned contrasts

#this is a bad function for bad people

max_var <- function(x,n){
  combinations <- t(combn(x,n))
  indexes <- t(combn(1:length(x),n))
  vars <- apply(combinations,1,var)
  return(indexes[which.max(vars),])
}
emm1 = emmeans(mod_nw_pre_36, specs = ~ S0_ID*Treatment_history)
emm1
emms <- summary(emm1)
x=emms$emmean[!is.na(emms$emmean)]
indexes=max_var(x,4)
emms$Treatment_history[indexes]

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
#write.csv(adjust.p.df, "Results/SLA_nw_planned.csv")
#nothing

vis_sla_nw<-visreg(mod_nw, xvar="Treatment_history", by="PrePeak") 
Res_sla_nw<-vis_sla_nw$res # Extract residuals

Res_sla_nw$PrePeak <- factor(Res_sla_nw$PrePeak, level = c('Pre', 'Peak'))
Res_sla_nw$Treatment_history <- factor(Res_sla_nw$Treatment_history, level = c('DD', 'DW', "WD", "WW"), labels = c("DDW", "DWW", "WDW", "WWW"))


plot_nw <- ggplot(Res_sla_nw, aes(Treatment_history, y=visregRes))+
  geom_violin(aes(), fill ="skyblue3", position = position_dodge(width = 1))+
  stat_summary(aes(), fun=mean, colour="black",position = position_dodge(width = 1))+
  stat_summary(aes(),fun.data=mean_se, geom="errorbar", colour="black", width=0.5,position = position_dodge(width = 1))+
  xlab("Treatment History") +
  ylim(-30,700)+
  facet_wrap(.~ PrePeak) +
  ylab("Specific Leaf Area")+
  scale_color_manual(values= c("Dry"="#FFA100", "Wet"="skyblue3"))+
  scale_fill_manual(values= c("Dry"="#FFA100", "Wet"="skyblue3"))+
  theme_classic()
plot_nw <-plot_nw + theme(
  axis.text.x = element_text(size=12, face="bold", angle=0,hjust=0.5),
  axis.text.y = element_text(size=15,face="bold"),
  axis.title.x = element_text(color="black", size=20, vjust = 0.5, face="bold"),
  axis.title.y = element_text(color="black", size=20,vjust = 2, face="bold",hjust=0.5))
plot_nw <-plot_nw + facet_wrap(.~ PrePeak) +
  theme(legend.title = element_blank(),legend.text = element_text(size=12,face="bold"),
        strip.background = element_blank(), strip.text.x=element_text(size=14,face="bold",hjust=0.05,vjust=-1.2))
plot_nw


###### DRY ############

## 
mod_nd <- lmer(sla~Rep + Site + PrePeak*Treatment_history + (1|S0_ID), data=dat_north_dry)
#singularity 

#simpler model 
mod_nd2 <-  lm(sla~Rep + Site + PrePeak*Treatment_history, data=dat_north_dry)
#models show same results, so singularity alert is okay

#check for model violations
plot(mod_nd, type=c("p","smooth"), col.line=1) #fitted vs residuals

plot(mod_nd,
     sqrt(abs(resid(.)))~fitted(.),
     type=c("p","smooth"), col.line=1) #scale-location plot

lattice::qqmath(mod_nd) # qqplot
plot(mod_nd2)

hist(resid(mod_nd2)) #histogram

plot(mod_nd, rstudent(.) ~ hatvalues(.)) #residuals vs leverage

#anova 
anova_nd <- anova(mod_nd)
anova_nd
# Treatment history p = 0.01
# Year*TrtHist p = 0.03
#write.csv(anova_nd, "Results/sla_nd_anova.csv")


#grab effect sizes of the fixed factors
#es = difference between groups/standard deviation of the variance 
effsize_nd <- omega_squared(mod_nd)
effsize_nd
#TRTHIS = 0.02
#PrePeak*TRTHIS = 0.02
#write.csv(effsize_nd, "Results/sla_nd_effectsize.csv")


###
emm1 = emmeans(mod_nd, specs = ~ PrePeak*Treatment_history)
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
#write.csv(adjust.p.df, "Results/sla_nd_planned.csv")
#peak DD-WW = 0.014 --> 
#peak DD-DW = 0.03 -->
#pre DW-WW = 0.05 -->
#peak DD-WD = 0.08 -->
#pre-peak DD = 0.09 -->
#pre-peak DW = 0.09 -->

vis_sla_nd<-visreg(mod_nd, xvar="Treatment_history", by="PrePeak") 
Res_sla_nd<-vis_sla_nd$res # Extract residuals

Res_sla_nd$PrePeak <- factor(Res_sla_nd$PrePeak, level = c('Pre', 'Peak'))
Res_sla_nd$Treatment_history <- factor(Res_sla_nd$Treatment_history, level = c('DD', 'DW', "WD", "WW"), labels = c("DDD", "DWD", "WDD", "WWD"))


plot_nd <- ggplot(Res_sla_nd, aes(Treatment_history, y=visregRes))+
  geom_violin(aes(), fill ="#FFA100", position = position_dodge(width = 1))+
  stat_summary(aes(), fun=mean, colour="black",position = position_dodge(width = 1))+
  stat_summary(aes(),fun.data=mean_se, geom="errorbar", colour="black", width=0.5,position = position_dodge(width = 1))+
  xlab("Treatment History") +
  ylim(-30,700)+
  facet_wrap(.~ PrePeak) +
  ylab("Specific Leaf Area")+
  scale_color_manual(values= c("Dry"="#FFA100", "Wet"="skyblue3"))+
  scale_fill_manual(values= c("Dry"="#FFA100", "Wet"="skyblue3"))+
  theme_classic()
plot_nd <-plot_nd + theme(
  axis.text.x = element_text(size=12, face="bold", angle=0,hjust=0.5),
  axis.text.y = element_text(size=15,face="bold"),
  axis.title.x = element_text(color="black", size=20, vjust = 0.5, face="bold"),
  axis.title.y = element_text(color="black", size=20,vjust = 2, face="bold",hjust=0.5))
plot_nd <-plot_nd + facet_wrap(.~ PrePeak) +
  theme(legend.title = element_blank(),legend.text = element_text(size=12,face="bold"),
        strip.background = element_blank(), strip.text.x=element_text(size=14,face="bold",hjust=0.05,vjust=-1.2))
plot_nd
#save 8x6


################### SOUTH ############################################
#### WET ######
mod_sw <- lmer(sla~Rep + Site + PrePeak*Treatment_history + (1|S0_ID), data=dat_south_wet)

#check for model violations
plot(mod_sw, type=c("p","smooth"), col.line=1) #fitted vs residuals

plot(mod_sw,
     sqrt(abs(resid(.)))~fitted(.),
     type=c("p","smooth"), col.line=1) #scale-location plot

lattice::qqmath(mod_sw) # qqplot

hist(resid(mod_sw)) #histogram

plot(mod_sw, rstudent(.) ~ hatvalues(.)) #residuals vs leverage

#anova
anova_sw <- anova(mod_sw) 
anova_sw
#no differences 
#write.csv(anova_sw, "Results/sla_sw_anova.csv")


#grab effect sizes of the fixed factors
#es = difference between groups/standard deviation of the variance 
effsize_sw <- omega_squared(mod_sw)
effsize_sw
#write.csv(effsize_sw, "Results/sla_sw_effectsize.csv")

###

emm1 = emmeans(mod_sw, specs = ~ PrePeak*Treatment_history)
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
#no differences
#write.csv(adjust.p.df, "Results/sla_sw_planned.csv")

vis_sla_sw<-visreg(mod_sw, xvar="Treatment_history", by="PrePeak") 
Res_sla_sw<-vis_sla_sw$res # Extract residuals

Res_sla_sw$PrePeak <- factor(Res_sla_sw$PrePeak, level = c('Pre', 'Peak'))
Res_sla_sw$Treatment_history <- factor(Res_sla_sw$Treatment_history, level = c('DD', 'DW', "WD", "WW"), labels = c("DDW", "DWW", "WDW", "WWW"))


plot_sw <- ggplot(Res_sla_sw, aes(Treatment_history, y=visregRes))+
  geom_violin(aes(), fill ="skyblue3", position = position_dodge(width = 1))+
  stat_summary(aes(), fun=mean, colour="black",position = position_dodge(width = 1))+
  stat_summary(aes(),fun.data=mean_se, geom="errorbar", colour="black", width=0.5,position = position_dodge(width = 1))+
  xlab("Treatment History") +
  ylim(-30,700)+
  facet_wrap(.~ PrePeak) +
  ylab("Specific Leaf Area")+
  scale_color_manual(values= c("Dry"="#FFA100", "Wet"="skyblue3"))+
  scale_fill_manual(values= c("Dry"="#FFA100", "Wet"="skyblue3"))+
  theme_classic()
plot_sw <-plot_sw + theme(
  axis.text.x = element_text(size=12, face="bold", angle=0,hjust=0.5),
  axis.text.y = element_text(size=15,face="bold"),
  axis.title.x = element_text(color="black", size=20, vjust = 0.5, face="bold"),
  axis.title.y = element_text(color="black", size=20,vjust = 2, face="bold",hjust=0.5))
plot_sw <-plot_sw + facet_wrap(.~ PrePeak) +
  theme(legend.title = element_blank(),legend.text = element_text(size=12,face="bold"),
        strip.background = element_blank(), strip.text.x=element_text(size=14,face="bold",hjust=0.05,vjust=-1.2))
plot_sw


###### DRY ############

## 
mod_sd <- lmer(sla~Rep + Site + PrePeak*Treatment_history + (1|S0_ID), data=dat_south_dry)

#check for model violations
plot(mod_sd, type=c("p","smooth"), col.line=1) #fitted vs residuals

plot(mod_sd,
     sqrt(abs(resid(.)))~fitted(.),
     type=c("p","smooth"), col.line=1) #scale-location plot

lattice::qqmath(mod_sd) # qqplot

hist(resid(mod_sd)) #histogram

plot(mod_sd, rstudent(.) ~ hatvalues(.)) #residuals vs leverage

#anova
aov_sla_sd <- anova(mod_sd)
aov_sla_sd
#write.csv(aov_sla_sd, "Results/sla_sd_anova.csv")
# Year p = 0.04

#grab effect sizes of the fixed factors
#es = difference between groups/standard deviation of the variance 
effsize_sd <- omega_squared(mod_sd)
effsize_sd
#PrePeak = 0.19
#write.csv(effsize_sd, "Results/sla_sd_effectsize.csv")

###
emm1 = emmeans(mod_sd, specs = ~ PrePeak*Treatment_history)
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
#write.csv(adjust.p.df, "Results/sla_sd_planned.csv")
#nothing

vis_sla_sd<-visreg(mod_sd, xvar="Treatment_history", by="PrePeak") 
Res_sla_sd<-vis_sla_sd$res # Extract residuals

Res_sla_sd$PrePeak <- factor(Res_sla_sd$PrePeak, level = c('Pre', 'Peak'))
Res_sla_sd$Treatment_history <- factor(Res_sla_sd$Treatment_history, level = c('DD', 'DW', "WD", "WW"), labels = c("DDD", "DWD", "WDD", "WWD"))


plot_sd <- ggplot(Res_sla_sd, aes(Treatment_history, y=visregRes))+
  geom_violin(aes(), fill ="#FFA100", position = position_dodge(width = 1))+
  stat_summary(aes(), fun=mean, colour="black",position = position_dodge(width = 1))+
  stat_summary(aes(),fun.data=mean_se, geom="errorbar", colour="black", width=0.5,position = position_dodge(width = 1))+
  xlab("Treatment History") +
  ylim(-30,700)+
  facet_wrap(.~ PrePeak) +
  ylab("Specific Leaf Area")+
  scale_color_manual(values= c("Dry"="#FFA100", "Wet"="skyblue3"))+
  scale_fill_manual(values= c("Dry"="#FFA100", "Wet"="skyblue3"))+
  theme_classic()
plot_sd <-plot_sd + theme(
  axis.text.x = element_text(size=12, face="bold", angle=0,hjust=0.5),
  axis.text.y = element_text(size=15,face="bold"),
  axis.title.x = element_text(color="black", size=20, vjust = 0.5, face="bold"),
  axis.title.y = element_text(color="black", size=20,vjust = 2, face="bold",hjust=0.5))
plot_sd <-plot_sd + facet_wrap(.~ PrePeak) +
  theme(legend.title = element_blank(),legend.text = element_text(size=12,face="bold"),
        strip.background = element_blank(), strip.text.x=element_text(size=14,face="bold",hjust=0.05,vjust=-1.2))
plot_sd
