#################### Biomass ##################################################

library(tidyverse)
library(lme4)
library(ggplot2)
library(visreg)
library(emmeans)
library(effectsize)
library(multtest)

rep1 <- read.csv("Data/S2_root_rep1.csv")
#rep1 <- subset(rep1, select=-c(5:11))
#rep1$Rep <- "1"
rep2 <- read.csv("Data/S2_root_rep2.csv")
#rep2 <- subset(rep2, select=-c(5:7))
#rep2$Rep <- "2"
#rep2$ww <- NULL 
dat <- rbind(rep1,rep2)

s0 <- read.csv("Data/s2_s0_ids.csv")
trthis <- read.csv("Data/Trthis_IDs.csv")

### merge dataframes together
dat2 <- merge(s0,trthis, by="S0_ID")
dat2 <- merge(dat,dat2, by="ID")

## merge S0 and S1 treatments
dat2 <- dat2 %>%
  rowwise()%>%
  mutate(Treatment_history=paste(S0_trt,S1_trt,sep=""))

dat2$dw<- as.numeric(dat2$dw)

dat2$Rep <- as.character(dat2$Rep)

#################### root biomass ##################################################

# break up the data by region and treatment 
dat_north_wet <- filter(dat2, Region == "North", Treatment == "Wet")
dat_north_dry <- filter(dat2, Region == "North", Treatment == "Dry")
dat_south_wet <- filter(dat2, Region == "South", Treatment == "Wet")
dat_south_dry <- filter(dat2, Region == "South", Treatment == "Dry")



################### NORTH ############################################

###### DRY ############

## 
mod_nd <- lmer(dw~Rep + Site + PrePeak*Treatment_history + (1|S0_ID), data=dat_north_dry)

dat_north_dry$dw_log <- log(dat_north_dry$dw)
mod_nd_log <- lmer(dw_log~Rep + Site + PrePeak*Treatment_history + (1|S0_ID), data=dat_north_dry)


#check for model violations
plot(mod_nd, type=c("p","smooth"), col.line=1) #fitted vs residuals

plot(mod_nd,
     sqrt(abs(resid(.)))~fitted(.),
     type=c("p","smooth"), col.line=1) #scale-location plot

lattice::qqmath(mod_nd) # qqplot
lattice::qqmath(mod_nd_log) # qqplot

hist(resid(mod_nd_log)) #histogram

plot(mod_nd, rstudent(.) ~ hatvalues(.)) #residuals vs leverage

#anova
anova_nd<- anova(mod_nd_log) 
anova_nd
# no sign.
#write.csv(anova_nd, "Results/root_nd_anova.csv")


###
emm1 = emmeans(mod_nd_log, specs = ~ PrePeak*Treatment_history)
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
#write.csv(adjust.p.df, "Results/root_nd_planned.csv")
#nothing

vis_fl_nd<-visreg(mod_nd_log, xvar="Treatment_history", by="PrePeak") 
Res_fl_nd<-vis_fl_nd$res # Extract residuals

Res_fl_nd$PrePeak <- factor(Res_fl_nd$PrePeak, level = c('Pre', 'Peak'))
Res_fl_nd$Treatment_history <- factor(Res_fl_nd$Treatment_history, level = c('DD', 'DW', "WD", "WW"), labels = c("DDD", "DWD", "WDD", "WWD"))


plot_nd <- ggplot(Res_fl_nd, aes(Treatment_history, y=exp(visregRes)))+
  geom_violin(aes(), fill ="#FFA100", position = position_dodge(width = 1))+
  stat_summary(aes(), fun=mean, colour="black",position = position_dodge(width = 1))+
  stat_summary(aes(),fun.data=mean_se, geom="errorbar", colour="black", width=0.5,position = position_dodge(width = 1))+
  xlab("Treatment History") +
 # facet_wrap(.~ PrePeak) +
  scale_y_continuous(name="Root Biomass (g)")+
  scale_color_manual(values= c("Dry"="#FFA100", "Wet"="skyblue3"))+
  scale_fill_manual(values= c("Dry"="#FFA100", "Wet"="skyblue3"))+
  theme_classic()
plot_nd <-plot_nd + theme(
  axis.text.x = element_text(size=12, face="bold", angle=0,hjust=0.5),
  axis.text.y = element_text(size=15,face="bold"),
  axis.title.x = element_text(color="black", size=20, vjust = 0.5, face="bold"),
  axis.title.y = element_text(color="black", size=20,vjust = 2, face="bold",hjust=0.5))
plot_nd <-plot_nd + #facet_wrap(.~ PrePeak) +
  theme(legend.title = element_blank(),legend.text = element_text(size=12,face="bold"),
        strip.background = element_blank(), strip.text.x=element_text(size=14,face="bold",hjust=0.05,vjust=-1.2))
plot_nd


################### SOUTH ############################################

###### DRY ############

## 
mod_sd <- lmer(dw~Rep + Site + PrePeak*Treatment_history + (1|S0_ID), data=dat_south_dry)

dat_south_dry$dw_log <- log(dat_south_dry$dw)
mod_sd_log <- lmer(dw_log~Rep + Site + PrePeak*Treatment_history + (1|S0_ID), data=dat_south_dry)

#check for model violations
plot(mod_sd, type=c("p","smooth"), col.line=1) #fitted vs residuals

plot(mod_sd,
     sqrt(abs(resid(.)))~fitted(.),
     type=c("p","smooth"), col.line=1) #scale-location plot

lattice::qqmath(mod_sd) # qqplot
lattice::qqmath(mod_sd_log) # qqplot

hist(resid(mod_sd)) #histogram
hist(resid(mod_sd_log)) #histogram

plot(mod_sd, rstudent(.) ~ hatvalues(.)) #residuals vs leverage

#anova
anova_sd <- anova(mod_sd_log) 
anova_sd
#no significance 
#write.csv(anova_sd, "Results/root_sd_anova.csv")


###
emm1 = emmeans(mod_sd_log, specs = ~ PrePeak*Treatment_history)
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
#write.csv(adjust.p.df, "Results/root_sd_planned.csv")
#nothing

vis_fl_sd<-visreg(mod_sd_log, xvar="Treatment_history", by="PrePeak") 
Res_fl_sd<-vis_fl_sd$res # Extract residuals

Res_fl_sd$PrePeak <- factor(Res_fl_sd$PrePeak, level = c('Pre', 'Peak'))
Res_fl_sd$Treatment_history <- factor(Res_fl_sd$Treatment_history, level = c('DD', 'DW', "WD", "WW"), labels = c("DDD", "DWD", "WDD", "WWD"))


plot_sd <- ggplot(Res_fl_sd, aes(Treatment_history, y=exp(visregRes)))+
  geom_violin(aes(), fill ="#FFA100", position = position_dodge(width = 1))+
  stat_summary(aes(), fun=mean, colour="black",position = position_dodge(width = 1))+
  stat_summary(aes(),fun.data=mean_se, geom="errorbar", colour="black", width=0.5,position = position_dodge(width = 1))+
  xlab("Treatment History") +
 # facet_wrap(.~ PrePeak) +
  scale_y_continuous(name="Root Biomass (g)")+
  scale_color_manual(values= c("Dry"="#FFA100", "Wet"="skyblue3"))+
  scale_fill_manual(values= c("Dry"="#FFA100", "Wet"="skyblue3"))+
  theme_classic()
plot_sd <-plot_sd + theme(
  axis.text.x = element_text(size=12, face="bold", angle=0,hjust=0.5),
  axis.text.y = element_text(size=15,face="bold"),
  axis.title.x = element_text(color="black", size=20, vjust = 0.5, face="bold"),
  axis.title.y = element_text(color="black", size=20,vjust = 2, face="bold",hjust=0.5))
plot_sd <-plot_sd + # facet_wrap(.~ PrePeak) +
  theme(legend.title = element_blank(),legend.text = element_text(size=12,face="bold"),
        strip.background = element_blank(), strip.text.x=element_text(size=14,face="bold",hjust=0.05,vjust=-1.2))
plot_sd
